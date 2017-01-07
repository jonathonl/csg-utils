#!/bin/bash
set -euo pipefail


for sample_id in "NWD000"
do
  zones=("us-central1-a" "us-central1-b" "us-central1-c" "us-central1-f")
  machine_zone=${zones[$[ $RANDOM % ${#zones[@]} ]]}

  machine_name=$(echo "pre-align-"$(basename $sample_id | cut -f 1 -d '.') | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
  #MACHINE_TYPE_OPTS="--machine-type n1-standard-2"
  machine_type_opts="--custom-cpu 1 --custom-memory 6656MiB"

  gcloud compute instances create --zone $machine_zone $machine_name --scopes cloud-platform --image ubuntu-1604-alignment $machine_type_opts --boot-disk-size 300 --metadata "sample-id=${sample_id}" --metadata-from-file "startup-script=pre-align.sh" &
  sleep 5s
done

wait
