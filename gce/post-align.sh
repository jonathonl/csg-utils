#!/bin/bash
set -uo pipefail

export HOME=/root

ulimit -n 8192

run()
{
  local sample_id=$1
  [[ -z $sample_id ]] && echo "sample_id is empty" && return -1

  echo "[$(date)] Downloading input cram files"
  start_time=$(date +%s)
  gsutil -q cp "gs://topmed-crams/${sample_id}/*.cram" /home/alignment/
  rc=$?
  echo "[$(date)] Download exit status: ${rc}"
  echo "[$(date)] Download elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  echo "[$(date)] Starting pipeline"
  start_time=$(date +%s)

  input_crams_read_count=0
  for input_file in /home/alignment/*.cram 
  do 
    echo "[$(date)] Sorting ${input_file} ..."
    samtools flagstat $input_file > ${input_file}.flagstat
    input_crams_read_count=$(( $input_crams_read_count + $(grep 'paired in sequencing' ${input_file}.flagstat | awk '{print $1}') ))
    tmp_prefix=${input_file%.cram}.tmp
    samtools sort --reference /home/alignment/ref/hs38DH.fa --threads 1 -T $tmp_prefix -o ${input_file%.cram}.sorted.bam $input_file
    rc=$?
    [[ $rc != 0 ]] && break
    rm -f $input_file ${input_file}.flagstat ${tmp_prefix}*
  done

  if [[ $rc == 0 ]]
  then 
    echo "[$(date)] Merging ..."
    samtools merge --threads 1 /home/alignment/merged.bam /home/alignment/*.sorted.bam \
    && rm /home/alignment/*.sorted.bam \
    && bam-non-primary-dedup dedup_LowMem --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30 --log /home/alignment/dedup_lowmem.metrics --recab --in /home/alignment/merged.bam --out -.ubam --refFile /home/alignment/ref/hs38DH.fa --dbsnp /home/alignment/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
      | samtools view -h -C -T /home/alignment/ref/hs38DH.fa -o /home/alignment/output.cram --threads 1
    rc=$?
    
    if [[ $rc == 0 ]]
    then
      unique_sample_count=$(samtools view -H /home/alignment/output.cram | grep -o "\sSM:NWD[0-9]*" | sort | uniq | wc -l)
      if [[ $unique_sample_count != 1 ]]
      then
        echo "[$(date)] Multiple samples found in output (${unique_sample_count})."
        rc=-1
      fi
      
      samtools flagstat /home/alignment/output.cram > /home/alignment/output.cram.flagstat
      output_cram_read_count=$(grep 'paired in sequencing' /home/alignment/output.cram.flagstat | awk '{print $1}')

      echo "[$(date)] Input read count: ${input_crams_read_count}"
      echo "[$(date)] Output read count: ${output_cram_read_count}"

      if [[ $input_crams_read_count != $output_cram_read_count || $output_cram_read_count == 0 ]] 
      then 
        echo "[$(date)] Failed flagstat validation."
        rc=-1
      fi
    fi
  fi
  echo "[$(date)] Pipeline exit status: ${rc}"
  echo "[$(date)] Pipeline elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  output_uri=gs://topmed-recabs/${sample_id}/${sample_id}.recab.cram
  echo "[$(date)] Uploading ouput cram (${output_uri})"
  start_time=$(date +%s)
  
  gsutil -q cp /home/alignment/output.cram $output_uri && gsutil -q cp /home/alignment/output.cram.flagstat $output_uri".flagstat"
  rc=$?
  echo "[$(date)] Upload exit status: ${rc}"
  echo "[$(date)] Upload elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  echo "[$(date)] Deleting input cram files"
  gsutil rm gs://topmed-crams/${sample_id}/*
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
  next_sample=$(mysql -NB topmed_remapping -e "SELECT samples.id FROM samples LEFT JOIN statuses ON samples.status_id=statuses.id WHERE samples.compute_node_id='${compute_node_id}' AND statuses.name='running-post-align' ORDER BY last_updated DESC LIMIT 1")

  while [[ $continue_running == 1 ]]
  do
    if [[ -z $next_sample ]]
    then
      for i in {1..5}
      do
        next_sample=$(mysql -NB topmed_remapping -e "\
          START TRANSACTION; \
          SET @aligned_status_id = (SELECT id FROM statuses WHERE name='aligned'); \
          SET @running_post_align_status_id = (SELECT id FROM statuses WHERE name='running-post-align'); \
          UPDATE samples \
          SET compute_node_id='${compute_node_id}', status_id=@running_post_align_status_id \
          WHERE status_id=@aligned_status_id AND compute_node_id IS NULL \
          ORDER BY last_updated ASC LIMIT 1; \
          SELECT id FROM samples WHERE compute_node_id='${compute_node_id}' AND status_id=@running_post_align_status_id ORDER BY last_updated DESC LIMIT 1; \
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
      run_status=$( [[ $? == 0 ]] && echo "post-aligned" || echo "failed-post-align" )
      run_elapsed_seconds=$(( $(date +%s) - $run_start_time ))

      mysql topmed_remapping -e "INSERT INTO job_executions (compute_node_id, sample_id, start_date, elapsed_seconds) VALUES ('${compute_node_id}', '${next_sample}', FROM_UNIXTIME(${run_start_time}), ${run_elapsed_seconds})"

      for i in {1..5}
      do
        mysql topmed_remapping -e "UPDATE samples SET compute_node_id=NULL, status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${next_sample}'" && break || sleep $(( $i * 5 ))s
      done

      gzip < /home/alignment/run.log > /home/alignment/run.log.gz
      gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${next_sample}/post_align_${run_start_time}.log.gz
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
