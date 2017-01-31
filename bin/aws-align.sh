#!/bin/bash
RG_LINE=$1
shift
OUTPUT_DIR=$1
shift
INPUT_FILES=$@

set -o pipefail

MACHINE_NAME="align-"$(head -c16 /dev/urandom | xxd -p | tr -d \\n)
MACHINE_TAG=$(basename $1 .fastq.gz | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)

docker-machine create --driver amazonec2 --amazonec2-tags $MACHINE_TAG --amazonec2-ami ami-cc0579db --amazonec2-instance-type c4.8xlarge --amazonec2-root-size 300 --amazonec2-request-spot-instance --amazonec2-spot-price 0.50 $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
  COMMAND="set -o pipefail; rm -f /home/alignment/output.cram /home/alignment/output.cram.ok; bwa mem -t 36 -K 100000000 -Y -p -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 36 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - && touch /home/alignment/output.cram.ok"
  
  CONTAINER_ID=$(docker-machine ssh $MACHINE_NAME sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c \""$COMMAND"\")

  EXIT_STATUS=$?
  if [[ $EXIT_STATUS == 0 ]]
  then  
    for f in $INPUT_FILES
    do
      echo "Uploading "$f" ..."
      docker-machine scp $f $MACHINE_NAME":/home/alignment/input.fastq.gz"
      EXIT_STATUS=$?

      RETRY_COUNTER=0
      while [[ $EXIT_STATUS == 0 && $RETRY_COUNTER -lt 5 ]]
      do  

        docker-machine ssh $MACHINE_NAME sudo docker start $CONTAINER_ID
        EXIT_STATUS=$?

        if [[ $EXIT_STATUS == 0 ]]
        then
          CONTAINER_IS_RUNNING=1
          MACHINE_IS_RUNNING=1
          FAILED_CONTAINER_POLL_COUNT=0
          while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
          do
            sleep 60s
            
            CONTAINER_STATUS=$(docker-machine ssh $MACHINE_NAME sudo docker ps -al --format {{.Status}})
            if [[ $CONTAINER_STATUS =~ Exited\ \((.*)\) ]]
            then
              CONTAINER_IS_RUNNING=0
              EXIT_STATUS=${BASH_REMATCH[1]}
              
              echo "Fetching logs ..."
              docker-machine ssh $MACHINE_NAME "sudo docker logs $CONTAINER_ID"
              if [[ $? != 0 ]]
              then
                echo "Fetching logs FAILED!"
              fi
            fi

            if [[ $CONTAINER_STATUS ]]
            then
              FAILED_CONTAINER_POLL_COUNT=0
            else
              let FAILED_CONTAINER_POLL_COUNT++
              echo "Failed Count: "$FAILED_CONTAINER_POLL_COUNT
            fi

          done

          if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
          then
            OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $f .fastq.gz)".cram"
            echo "Downloading "$OUTPUT_FILE" ..."
            docker-machine scp $MACHINE_NAME":/home/alignment/output.cram" $OUTPUT_FILE && docker-machine scp $MACHINE_NAME":/home/alignment/output.cram.ok" $OUTPUT_FILE".ok"
            EXIT_STATUS=$?
            break
          elif [[ $FAILED_CONTAINER_POLL_COUNT  == 5 ]] 
          then
            EXIT_STATUS=-1
            echo "Machine stopped"
            #sleep 300s
            #echo "Restarting "$MACHINE_NAME" ..."
            #docker-machine start $MACHINE_NAME
          fi
        fi

        let RETRY_COUNTER++
      done

      if [[ $EXIT_STATUS != 0 ]]
      then
        break
      fi

    done
  fi
fi

# Machine may exist even if "docker-machine create" fails.
docker-machine rm --force $MACHINE_NAME

exit $EXIT_STATUS
