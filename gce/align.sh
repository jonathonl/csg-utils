#!/bin/bash
set -uo pipefail

export HOME=/root

sample_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/sample-id" -H "Metadata-Flavor: Google" --silent)

run()
{
  [[ -z $sample_id ]] && echo "sample_id is empty" && return -1

  echo "[$(date)] Parsing todo"
  start_time=$(date +%s)
  
  gsutil ls gs://topmed-crams/${sample_id}/*.cram.ok | xargs -I % basename % .cram.ok > /home/alignment/${sample_id}.completed-list
  gsutil -q cp gs://topmed-fastqs/${sample_id}/${sample_id}.list /home/alignment/${sample_id}.list \
  && grep -v -f /home/alignment/${sample_id}.completed-list /home/alignment/${sample_id}.list > /home/alignment/${sample_id}.todo-list \
  && todo_count=$(( $(cat /home/alignment/${sample_id}.todo-list | wc -l) - 1 ))
  rc=$?
  echo "[$(date)] Parsing todo exit status: ${rc}"
  echo "[$(date)] Parsing todo elapsed time: "$(( $(date +%s) - $start_time ))"s"
  
  [[ $rc != 0 ]] && return $rc

  echo "Todo count: "$todo_count  

  if [[ $todo_count -gt 0 ]]
  then
    while read line
    do
      line_rg=$(echo $line | cut -d ' ' -f 4- | sed -e "s/ /\\\t/g")
      
      input_file_name=$(basename $(echo $line | cut -f 2 -d ' '))
      if [[ -z $input_file_name ]]
      then
        echo "input_file_name is empty"
        rc=-1
        break
      fi
      
      input_uri=gs://topmed-fastqs/${sample_id}/$input_file_name

      #gce-align.sh 1 $sample_id $rg_counter "$line_rg" gs://topmed-crams/${sample_id} $input_files
      echo "[$(date)] Downloading input cram (${input_uri})"
      start_time=$(date +%s)
      gsutil -q cp $input_uri /home/alignment/input.fastq.gz
      rc=$?
      echo "[$(date)] Download exit status: ${rc}"
      echo "[$(date)] Download elapsed time: "$(( $(date +%s) - $start_time ))"s"

      if [[ $rc == 0 ]]
      then
        echo "[$(date)] Starting pipeline"
        start_time=$(date +%s)

        rm -f /home/alignment/output.cram /home/alignment/output.cram.ok \
        && bwa mem -t 32 -K 100000000 -Y -p -R "$line_rg" /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 32 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - \
        && touch /home/alignment/output.cram.ok

        rc=$?
        echo "[$(date)] Pipeline exit status: ${rc}"
        echo "[$(date)] Pipeline elapsed time: "$(( $(date +%s) - $start_time ))"s"

        if [[ $rc == 0 ]]
        then
          output_uri="gs://topmed-crams/${sample_id}/"$(basename $input_file_name .fastq.gz)".cram"
          echo "[$(date)] Uploading ouput cram (${output_uri})"
          start_time=$(date +%s)

          gsutil -q cp /home/alignment/output.cram $output_uri && gsutil -q cp /home/alignment/output.cram.ok $output_uri".ok"
          rc=$?
          echo "[$(date)] Upload exit status: ${rc}"
          echo "[$(date)] Upload elapsed time: "$(( $(date +%s) - $start_time ))"s"
        fi
      fi

      [[ $rc != 0 ]] && break
    done <<< "$(tail -n +2 /home/alignment/${sample_id}.todo-list)"
  fi

  [[ $rc != 0 ]] && return $rc


  echo "[$(date)] Deleting input fastq files"
  echo "gsutil rm gs://topmed-fastqs/${sample_id}/*"
  rc=$?
  echo "[$(date)] Delete exit status: ${rc}"

  return $rc
}

if [[ $sample_id ]]
then
  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='running-align') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  run &>> /home/alignment/run.log
  run_status=$( [[ $? == 0 ]] && echo "aligned" || echo "failed-align" )

  for i in {1..5}
  do
    mysql topmed_remapping -e "UPDATE samples SET status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${sample_id}'" && break || sleep $(( $i * 5 ))s
  done

  gzip /home/alignment/run.log
  gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${sample_id}/align_$(date +%s).log.gz
else
  echo "sample_id is empty"
fi

for i in {1..5}
do
  gcloud compute instances delete $(hostname) --quiet && break || sleep $(( $i * 5 ))s
done