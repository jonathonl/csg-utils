#!/bin/bash
set -uo pipefail

export HOME=/root

sample_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/sample-id" -H "Metadata-Flavor: Google" --silent)

run()
{
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
    | samtools sort -l 1 -@ 1 -m 3000M -n -T /home/alignment/sort_tmp - \
    | samtools fixmate - - \
    | bam-ext-mem-sort-manager bam2fastq --in -.bam --outBase $local_output_base --maxRecordLimitPerFq 20000000 --sortByReadNameOnTheFly --readname --gzip
  
  rc=$?
  echo "[$(date)] Pipeline exit status: ${rc}"
  echo "[$(date)] Pipeline elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  fastq_reads=$(( $(zcat /home/alignment/*.fastq.gz | wc -l) / 4 ));
  samtools flagstat $local_input_file > /home/alignment/cram_flagstat.txt;
  cram_reads=$(grep 'paired in sequencing' /home/alignment/cram_flagstat.txt | awk '{print $1}');
  
  echo '[$(date)] Cram read count: '$cram_reads;
  echo '[$(date)] Fastq read count: '$fastq_reads;

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

if [[ $sample_id ]]
then
  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='running-pre-align') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  run &> /home/alignment/run.log
  run_status=$( [[ $? == 0 ]] && echo "pre-aligned" || echo "failed-pre-align" )

  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  gzip /home/alignment/run.log
  gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${sample_id}/pre_align_$(date +%s).log.gz
else
  echo "sample_id is empty"
fi

for i in {1..5}
do
  gcloud compute instances delete $(hostname) --quiet && break || sleep $(( $i * 5 ))s
done