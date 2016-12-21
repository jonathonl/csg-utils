#!/bin/bash
#
# Usage:
#   gce-post-align.sh [input_uri] [output_uri]
#
#   [input_uri]  - Google Storage Bucket URI for aligned crams (e.g. gs://topmed-crams/<sample_id>
#   [output_uri] - Google Storage Bucket URI for recab crams (e.g. gs://topmed-recabs/<sample_id>/<recab_file>
#
set -u -o pipefail

ZONES=("us-central1-a" "us-central1-b" "us-central1-c" "us-central1-f")
MACHINE_ZONE=${ZONES[$[ $RANDOM % ${#ZONES[@]} ]]}

MACHINE_NAME=$(echo "post-align-"$(basename $1) | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
#MACHINE_TYPE_OPTS="--machine-type n1-standard-2"
MACHINE_TYPE_OPTS="--custom-cpu 1 --custom-memory 6656MiB"

gcloud compute instances create --zone $MACHINE_ZONE --scopes storage-full --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
  # Give ssh daemon some time to get running.
  sleep 30s

  SSH_COMMAND=$(gcloud compute ssh --dry-run --zone $MACHINE_ZONE $MACHINE_NAME)
  EXIT_STATUS=$?

  if [[ $EXIT_STATUS == 0 ]]
  then
    echo "[$(date)] Downloading input crams"
    START_TIME=$(date +%s)
    eval $SSH_COMMAND -- gsutil cp $1'/*.cram' /home/alignment/
    EXIT_STATUS=$?
    echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

    if [[ $EXIT_STATUS == 0 ]]
    then
      REMOTE_COMMAND="set -euo pipefail; for INPUT_FILE in /home/alignment/*.cram; do samtools sort --reference /home/alignment/ref/hs38DH.fa --threads 1 -T /home/alignment/tmp -o \\\${INPUT_FILE%.cram}.sorted.bam \\\$INPUT_FILE; rm \\\$INPUT_FILE; done; samtools merge --threads 1 /home/alignment/merged.bam /home/alignment/*.sorted.bam; rm /home/alignment/*.sorted.bam; bam-non-primary-dedup dedup_LowMem --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30 --log /home/alignment/dedup_lowmem.metrics --recab --in /home/alignment/merged.bam --out -.ubam --refFile /home/alignment/ref/hs38DH.fa --dbsnp /home/alignment/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz | samtools view -h -C -T /home/alignment/ref/hs38DH.fa -o /home/alignment/output.cram --threads 1; samtools flagstat /home/alignment/output.cram > /home/alignment/output.flagstat"

      echo "[$(date)] Creating container"
      CONTAINER_ID=$(eval $SSH_COMMAND -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c '\""$REMOTE_COMMAND"\"')

      echo "[$(date)] Starting container"
      eval $SSH_COMMAND -- sudo docker start $CONTAINER_ID

      CONTAINER_IS_RUNNING=1
      FAILED_CONTAINER_POLL_COUNT=0
      while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
      do
        sleep 180s

        CONTAINER_STATUS=$(eval $SSH_COMMAND -- sudo docker ps -al --format {{.Status}})
        if [[ $CONTAINER_STATUS =~ Exited\ \((.*)\) ]]
        then
          CONTAINER_IS_RUNNING=0
          EXIT_STATUS=${BASH_REMATCH[1]}

          echo "[$(date)] Fetching logs"
          eval $SSH_COMMAND -- sudo docker logs $CONTAINER_ID
          if [[ $? != 0 ]]
          then
            echo "[$(date)] Fetching logs FAILED!"
          fi
        fi

        if [[ $CONTAINER_STATUS ]]
        then
          FAILED_CONTAINER_POLL_COUNT=0
        else
          let FAILED_CONTAINER_POLL_COUNT++
          echo "[$(date)] Failed Count: $FAILED_CONTAINER_POLL_COUNT"
        fi
      done

      if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
      then
        OUTPUT_FILE=$2
        echo "[$(date)] Uploading $OUTPUT_FILE"
        START_TIME=$(date +%s)
        eval $SSH_COMMAND -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/output.cram $OUTPUT_FILE && eval $SSH_COMMAND -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/output.flagstat $OUTPUT_FILE".flagstat"
        EXIT_STATUS=$?
        echo "[$(date)] Upload exit status: $EXIT_STATUS"
        echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
      elif [[ $FAILED_CONTAINER_POLL_COUNT == 5 ]]
      then
        echo "[$(date)] Lost communication with $MACHINE_NAME"
        EXIT_STATUS=-1
      fi
    fi
  fi
fi

# Machine may exist even if "docker-machine create" fails.
echo "[$(date)] Cleaning up ..."
gcloud compute instances delete --zone $MACHINE_ZONE --quiet $MACHINE_NAME

exit $EXIT_STATUS
