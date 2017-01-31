#!/bin/bash
set -uo pipefail

export HOME=/root

#sample_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/sample-id" -H "Metadata-Flavor: Google" --silent)

run()
{
  local sample_id=$1
  [[ -z $sample_id ]] && echo "sample_id is empty" && return -1

  input_uri=$(gsutil ls gs://topmed-incoming/**/${sample_id}.src.cram) #$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/input-uri" -H "Metadata-Flavor: Google" --silent)

  [[ $? != 0 || -z $input_uri ]] && echo "input_uri is empty" && return -1

  ouput_uri=gs://topmed-fastqs/${sample_id}/ # !!! Trailing slash important.

  local_input_file="/home/alignment/"$(basename $input_uri)
  local_output_base="/home/alignment/"$(basename $input_uri | cut -f 1 -d '.')


  echo "[$(date)] Downloading input cram file (${input_uri})"
  start_time=$(date +%s)
  gsutil -q cp $input_uri $local_input_file
  rc=$?
  echo "[$(date)] Download exit status: ${rc}"
  echo "[$(date)] Download elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  echo "[$(date)] Starting pipeline"
  start_time=$(date +%s)
  export REF_CACHE=/home/alignment/ref/md5/%2s/%2s/%s &&
  samtools view -uh -F 0x900 $local_input_file \
    | bam-ext-mem-sort-manager squeeze --in -.ubam --keepDups --rmTags AS:i,BD:Z,BI:Z,XS:i,MC:Z,MD:Z,NM:i,MQ:i --out -.ubam \
    | samtools sort -l 1 -@ 1 -m 3000M -n -T /home/alignment/${sample_id}.samtools_sort_tmp - \
    | samtools fixmate - - \
    | bam-ext-mem-sort-manager bam2fastq --in -.bam --outBase $local_output_base --maxRecordLimitPerFq 20000000 --sortByReadNameOnTheFly --readname --gzip
  
  rc=$?
  echo "[$(date)] Pipeline exit status: ${rc}"
  echo "[$(date)] Pipeline elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  fastq_reads=$(( $(zcat /home/alignment/*.fastq.gz | wc -l) / 4 ))
  samtools flagstat $local_input_file > /home/alignment/cram_flagstat.txt
  cram_reads=$(grep 'paired in sequencing' /home/alignment/cram_flagstat.txt | awk '{print $1}')
  
  echo "[$(date)] Cram read count: ${cram_reads}"
  echo "[$(date)] Fastq read count: ${fastq_reads}"

  if [[ $fastq_reads != $cram_reads || $cram_reads == 0 ]] 
  then 
    echo "[$(date)] Failed flagstat validation."
    return -1
  fi

  echo "[$(date)] Uploading ouput fastq files"
  start_time=$(date +%s)
  
  gsutil -q cp ${local_output_base}*.fastq.gz $ouput_uri && gsutil -q cp ${local_output_base}.list $ouput_uri
  rc=$?
  echo "[$(date)] Upload exit status: ${rc}"
  echo "[$(date)] Upload elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  echo "[$(date)] Deleting input cram file (${input_uri})"
  gsutil rm $input_uri
  rc=$?
  echo "[$(date)] Delete exit status: ${rc}"

  return $rc
}

# There's a bug in gcloud where /etc/hostname doesn't get updated correctly.
compute_node_id=""
for i in {1..5}
do
  compute_node_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/hostname" -H "Metadata-Flavor: Google" --silent | cut -f 1 -d ".")
  [[ -z $compute_node_id ]] && sleep $(( $i * 5 ))s || break
done

if [[ -z $compute_node_id ]]
then
  echo "empty compute_node_id"
  exit -1
fi

continue_running=1

mysql topmed_remapping -e "INSERT INTO compute_nodes (id) VALUES ('${compute_node_id}') ON DUPLICATE KEY UPDATE enabled=IF(enabled != 0, 1, 0)"

if [[ $? == 0 ]]
then
  next_sample=$(mysql -NB topmed_remapping -e "SELECT samples.id FROM samples LEFT JOIN statuses ON samples.status_id=statuses.id WHERE samples.compute_node_id='${compute_node_id}' AND statuses.name='running-pre-align' ORDER BY last_updated DESC LIMIT 1")

  while [[ $continue_running == 1 ]]
  do
    if [[ -z $next_sample ]]
    then
      for i in {1..5}
      do
        next_sample=$(mysql -NB topmed_remapping -e "\
          START TRANSACTION; \
          SET @unprocessed_status_id = (SELECT id FROM statuses WHERE name='unprocessed'); \
          SET @running_pre_align_status_id = (SELECT id FROM statuses WHERE name='running-pre-align'); \
          UPDATE samples \
          SET compute_node_id='${compute_node_id}', status_id=@running_pre_align_status_id \
          WHERE status_id=@unprocessed_status_id AND compute_node_id IS NULL \
          ORDER BY last_updated ASC LIMIT 1; \
          SELECT id FROM samples WHERE compute_node_id='${compute_node_id}' AND status_id=@running_pre_align_status_id ORDER BY last_updated DESC LIMIT 1; \
          COMMIT;")
        [[ $? == 0 ]] && break || sleep $(( $i * 5 ))s
      done
    fi

    if [[ -z $next_sample ]]
    then
      echo "sample is empty"
      break
    else
      find /home/alignment -maxdepth 1 -type f | xargs rm -f
      run_start_time=$(date +%s)
      run $next_sample &> /home/alignment/run.log
      run_status=$( [[ $? == 0 ]] && echo "pre-aligned" || echo "failed-pre-align" )
      run_elapsed_seconds=$(( $(date +%s) - $run_start_time ))

      mysql topmed_remapping -e "INSERT INTO job_executions (compute_node_id, sample_id, start_date, elapsed_seconds) VALUES ('${compute_node_id}', '${next_sample}', FROM_UNIXTIME(${run_start_time}), ${run_elapsed_seconds})"

      for i in {1..5}
      do
        mysql topmed_remapping -e "UPDATE samples SET compute_node_id=NULL, status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${next_sample}'" && break || sleep $(( $i * 5 ))s
      done

      gzip < /home/alignment/run.log > /home/alignment/run.log.gz
      gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${next_sample}/pre_align_${run_start_time}.log.gz
    fi

    next_sample=""
    continue_running=$(mysql -NB topmed_remapping -e "SELECT enabled FROM compute_nodes WHERE id='${compute_node_id}'")
  done
fi

if [[ $continue_running == 2 ]]
then
  reboot
else
  for i in {1..10}
  do
    gcloud compute instances delete $compute_node_id --quiet && break || sleep $(( $i * 5 ))s
  done
fi
