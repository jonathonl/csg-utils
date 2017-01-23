#!/bin/bash
set -uo pipefail

export HOME=/root

#sample_id=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/sample-id" -H "Metadata-Flavor: Google" --silent)

run()
{
  local sample_id=$1
  [[ -z $sample_id ]] && echo "sample_id is empty" && return -1

  echo "[$(date)] Parsing todo"
  start_time=$(date +%s)
  
  rm -f /home/alignment/full_list.txt /home/alignment/completed_list.txt /home/alignment/todo_list.txt
  gsutil ls gs://topmed-crams/${sample_id}/*.cram.ok | xargs -I % basename % .cram.ok > /home/alignment/completed_list.txt
  gsutil -q cp gs://topmed-fastqs/${sample_id}/${sample_id}.list /home/alignment/full_list.txt \
  && grep -v -f /home/alignment/completed_list.txt /home/alignment/full_list.txt > /home/alignment/todo_list.txt \
  && todo_count=$(( $(cat /home/alignment/todo_list.txt | wc -l) - 1 ))
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
      echo "[$(date)] Downloading input fastq (${input_uri})"
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
          echo "[$(date)] Starting validation"
          start_time=$(date +%s)
          fastq_reads=$(( $(zcat /home/alignment/input.fastq.gz | wc -l) / 4 ))
          samtools flagstat /home/alignment/output.cram > /home/alignment/output.cram.flagstat
          cram_reads=$(grep 'paired in sequencing' /home/alignment/output.cram.flagstat | awk '{print $1}')
          
          echo "[$(date)] Cram read count: ${cram_reads}"
          echo "[$(date)] Fastq read count: ${fastq_reads}"

          if [[ $fastq_reads != $cram_reads || $cram_reads == 0 ]] 
          then 
            echo "[$(date)] Failed flagstat validation."
            rc=-1
          fi
          echo "[$(date)] Validation exit status: ${rc}"
          echo "[$(date)] Validation elapsed time: "$(( $(date +%s) - $start_time ))"s"


          if [[ $rc == 0 ]]
          then
            output_uri="gs://topmed-crams/${sample_id}/"$(basename $input_file_name .fastq.gz)".cram"
            echo "[$(date)] Uploading output cram (${output_uri})"
            start_time=$(date +%s)

            gsutil -q cp /home/alignment/output.cram $output_uri && gsutil -q cp /home/alignment/output.cram.ok $output_uri".ok"
            rc=$?
            echo "[$(date)] Upload exit status: ${rc}"
            echo "[$(date)] Upload elapsed time: "$(( $(date +%s) - $start_time ))"s"
          fi
        fi
      fi

      [[ $rc != 0 ]] && break
    done <<< "$(tail -n +2 /home/alignment/todo_list.txt)"
  fi

  [[ $rc != 0 ]] && return $rc


  echo "[$(date)] Deleting input fastq files"
  gsutil rm gs://topmed-fastqs/${sample_id}/*
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
  next_sample=$(mysql -NB topmed_remapping -e "SELECT samples.id FROM samples LEFT JOIN statuses ON samples.status_id=statuses.id WHERE samples.compute_node_id='${compute_node_id}' AND statuses.name='running-align' ORDER BY last_updated DESC LIMIT 1")

  while [[ $continue_running == 1 ]]
  do
    if [[ -z $next_sample ]]
    then
      for i in {1..5}
      do
        next_sample=$(mysql -NB topmed_remapping -e "\
          START TRANSACTION; \
          SET @pre_aligned_status_id = (SELECT id FROM statuses WHERE name='pre-aligned'); \
          SET @running_align_status_id = (SELECT id FROM statuses WHERE name='running-align'); \
          UPDATE samples \
          SET compute_node_id='${compute_node_id}', status_id=@running_align_status_id \
          WHERE status_id=@pre_aligned_status_id AND compute_node_id IS NULL \
          ORDER BY RAND() LIMIT 1; \
          SELECT id FROM samples WHERE compute_node_id='${compute_node_id}' AND status_id=@running_align_status_id ORDER BY last_updated DESC LIMIT 1; \
          COMMIT;")
        [[ $? == 0 ]] && break || sleep $(( $i * 5 ))s
      done

      # Just grabbed new sample, so start new log.
      find /home/alignment -maxdepth 1 -type f | xargs rm -f
    fi

    if [[ -z $next_sample ]]
    then
      echo "sample is empty"
      break
    else
      run_start_time=$(date +%s)
      run $next_sample &>> /home/alignment/run.log
      run_status=$( [[ $? == 0 ]] && echo "aligned" || echo "failed-align" )
      run_elapsed_seconds=$(( $(date +%s) - $run_start_time ))

      mysql topmed_remapping -e "INSERT INTO job_executions (compute_node_id, sample_id, start_date, elapsed_seconds) VALUES ('${compute_node_id}', '${next_sample}', FROM_UNIXTIME(${run_start_time}), ${run_elapsed_seconds})"

      for i in {1..5}
      do
        mysql topmed_remapping -e "UPDATE samples SET compute_node_id=NULL, status_id=(SELECT id FROM statuses WHERE name='${run_status}') WHERE id='${next_sample}'" && break || sleep $(( $i * 5 ))s
      done

      gzip < /home/alignment/run.log > /home/alignment/run.log.gz
      gsutil -q cp /home/alignment/run.log.gz gs://topmed-logs/${next_sample}/align_${run_start_time}.log.gz
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