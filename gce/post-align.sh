#!/bin/bash
set -uo pipefail

export HOME=/root

sample_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/sample-id" -H "Metadata-Flavor: Google" --silent)

run()
{
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

  for input_file in /home/alignment/*.cram 
  do 
    samtools sort --reference /home/alignment/ref/hs38DH.fa --threads 1 -T /home/alignment/tmp -o ${input_file%.cram}.sorted.bam $input_file
    rc=$?
    [[ $rc != 0 ]] && break
    rm $input_file
  done

  if [[ $rc == 0 ]]
  then 
    samtools merge --threads 1 /home/alignment/merged.bam /home/alignment/*.sorted.bam \
    && rm /home/alignment/*.sorted.bam \
    && bam-non-primary-dedup dedup_LowMem --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30 --log /home/alignment/dedup_lowmem.metrics --recab --in /home/alignment/merged.bam --out -.ubam --refFile /home/alignment/ref/hs38DH.fa --dbsnp /home/alignment/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
      | samtools view -h -C -T /home/alignment/ref/hs38DH.fa -o /home/alignment/output.cram --threads 1
    rc=$?
    
    if [[ $rc == 0 ]]
    then
      samtools flagstat /home/alignment/output.cram > /home/alignment/output.flagstat
      rc=$?
    fi
  fi
  echo "[$(date)] Pipeline exit status: ${rc}"
  echo "[$(date)] Pipeline elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  output_uri=gs://topmed-recabs/${sample_id}/${sample_id}.recab.cram
  echo "[$(date)] Uploading ouput cram (${output_uri})"
  start_time=$(date +%s)
  
  gsutil -q cp /home/alignment/output.cram $output_uri && gsutil -q cp /home/alignment/output.flagstat $output_uri".flagstat"
  rc=$?
  echo "[$(date)] Upload exit status: ${rc}"
  echo "[$(date)] Upload elapsed time: "$(( $(date +%s) - $start_time ))"s"

  [[ $rc != 0 ]] && return $rc

  echo "[$(date)] Deleting input cram files"
  echo "gsutil rm gs://topmed-crams/${sample_id}/*"
  rc=$?
  echo "[$(date)] Delete exit status: ${rc}"

  return $rc
}

if [[ $sample_id ]]
then
  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='running-post-align') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  run &> /home/alignment/run.log
  run_status=$( [[ $? == 0 ]] && echo "post-aligned" || echo "failed-post-align" )

  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  gzip /home/alignment/run.log
  gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${sample_id}/post_align_$(date +%s).log.gz
else
  echo "sample_id is empty"
fi

for i in {1..5}
do
  gcloud compute instances delete $(hostname) --quiet && break || sleep $(( $i * 5 ))s
done