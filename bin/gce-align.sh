#!/bin/bash
set -uo pipefail

#
# Usage:
#   gce-align.sh [paired] [sample_id] [read_group_ident] [read_group] [output_dir] [input_files]
#
#   [paired]           - either 0 or 1 for unpaired or paired reads
#   [sample_id]        - sample identifier (e.g. nwd123456)
#   [read_group_ident] - identification string to describe the read group for vm creatation (e.g rg-0)
#   [read_group]       - read group string from the .list file during bam2fastq
#   [output_dir]       - directory, or google bucket, where cram files will be written
#   [input_files]      - full path to fastq file(s). can be absolute path or google bucket uri (e.g. gs://topmed-fastqs/)
#

PAIRED=$1
shift
SAMPLE_ID=$1
shift
RG_IDENT=$1
shift
RG_LINE=$1
shift
OUTPUT_DIR=$1
shift
INPUT_FILES=$@
INPUT_FILES_COUNT=$#

case "$PAIRED" in
  0) BWA_PAIRED_OPT="" ;;
  1) BWA_PAIRED_OPT="-p" ;;
esac

ZONES=("us-central1-b" "us-central1-c" "us-central1-f")
MACHINE_ZONE=${ZONES[$[ $RANDOM % ${#ZONES[@]} ]]}

MACHINE_NAME=$(echo "align-${SAMPLE_ID}-${RG_IDENT}" | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
MACHINE_TYPE_OPTS="--custom-cpu 32 --custom-memory 64GiB"

EXIT_STATUS=0
gcloud compute instances create --zone $MACHINE_ZONE --scopes storage-full --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 --preemptible $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
# Give ssh daemon some time to get running.
  sleep 30s

  SSH_COMMAND=$(gcloud compute ssh --dry-run --zone $MACHINE_ZONE $MACHINE_NAME)
  EXIT_STATUS=$?

  if [[ $EXIT_STATUS == 0 ]]
  then
    COMMAND="set -o pipefail; rm -f /home/alignment/output.cram /home/alignment/output.cram.ok; bwa mem -t 32 -K 100000000 -Y $BWA_PAIRED_OPT -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 32 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - && touch /home/alignment/output.cram.ok"

    CONTAINER_ID=$(eval $SSH_COMMAND -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c '\""$COMMAND"\"')

    EXIT_STATUS=$?

    while [[ $# -gt 0 && $EXIT_STATUS == 0 ]]
    do
      F=$1
      echo "[$(date)] Downloading "$F" ..."
      START_TIME=$(date +%s)
      eval $SSH_COMMAND -- gsutil cp $F /home/alignment/input.fastq.gz
      EXIT_STATUS=$?
      echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

      if [[ $EXIT_STATUS == 0 ]]
      then
        eval $SSH_COMMAND -- sudo docker start $CONTAINER_ID
        EXIT_STATUS=$?

        if [[ $EXIT_STATUS == 0 ]] 
        then
          CONTAINER_IS_RUNNING=1
          FAILED_CONTAINER_POLL_COUNT=0
          while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
          do
            sleep 180s

            CONTAINER_STATUS_LINE=$(eval $SSH_COMMAND -- sudo docker ps -al --format {{.Status}})
            if [[ $CONTAINER_STATUS_LINE =~ Exited\ \((.*)\) ]]
            then
              CONTAINER_IS_RUNNING=0
              EXIT_STATUS=${BASH_REMATCH[1]}

              echo "[$(date)] Fetching logs ..."
              eval $SSH_COMMAND -- sudo docker logs $CONTAINER_ID
              if [[ $? != 0 ]]
              then
                echo "[$(date)] Fetching logs FAILED!"
              fi

              if [[ $EXIT_STATUS == 0 ]]
              then
                OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $F .fastq.gz)".cram"
                echo "[$(date)] Uploading "$OUTPUT_FILE" ..."
                START_TIME=$(date +%s)
                eval $SSH_COMMAND -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/output.cram $OUTPUT_FILE && eval $SSH_COMMAND -- gsutil cp /home/alignment/output.cram.ok $OUTPUT_FILE".ok"
                EXIT_STATUS=$?
                echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
                
                if [[ $EXIT_STATUS == 0 ]]
                then
                  shift #pops file off of input list
                fi
              fi
            fi

            if [[ $CONTAINER_STATUS_LINE ]]
            then
              FAILED_CONTAINER_POLL_COUNT=0
            else
              let FAILED_CONTAINER_POLL_COUNT++
              echo "[$(date)] Failed Count: "$FAILED_CONTAINER_POLL_COUNT
            fi
          done

          if [[ $FAILED_CONTAINER_POLL_COUNT == 5 ]]
          then
            echo "[$(date)] Machine stopped: "$MACHINE_NAME
            EXIT_STATUS=-1
          fi
        fi
      fi
    done
  fi
fi

# Machine may exist even if "docker-machine create" fails.
echo "[$(date)] Cleaning up ..."
gcloud compute instances delete --zone $MACHINE_ZONE --quiet $MACHINE_NAME


if [[ $# -gt 0 ]]
then
  echo "[$(date)] Giving up at "$1
else
  echo "[$(date)] Done."
fi

exit $EXIT_STATUS














