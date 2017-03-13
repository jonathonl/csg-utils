#!/bin/bash
set -euo pipefail

[[ $# -lt 1 ]] && exit -1

for node_int_id in "$@"
do
  zones=("us-central1-b" "us-central1-c" "us-central1-f")
  machine_zone=${zones[$[ $RANDOM % ${#zones[@]} ]]}

  machine_name="align-"$node_int_id
  machine_type_opts="--custom-cpu 32 --custom-memory 64GiB --preemptible"
  
  gcloud compute instances create --zone $machine_zone $machine_name --scopes cloud-platform --image ubuntu-1604-alignment $machine_type_opts --boot-disk-size 200 --metadata-from-file "startup-script=./align.sh" &
  sleep 5s
done

wait