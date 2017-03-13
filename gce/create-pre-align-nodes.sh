#!/bin/bash
set -euo pipefail

[[ $# -lt 1 ]] && exit -1

for node_int_id in "$@"
do
  zones=("us-central1-a" "us-central1-b" "us-central1-c" "us-central1-f")
  machine_zone=${zones[$[ $RANDOM % ${#zones[@]} ]]}

  machine_name="pre-align-"$node_int_id
  machine_type_opts="--custom-cpu 1 --custom-memory 6656MiB"
  
  gcloud compute instances create --zone $machine_zone $machine_name --scopes cloud-platform --image ubuntu-1604-alignment $machine_type_opts --boot-disk-size 300 --metadata-from-file "startup-script=./pre-align.sh" &
  sleep 5s
done

wait