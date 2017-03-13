#!/bin/bash

# example crontab
#SHELL=/bin/bash
#PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games
#
## m h  dom mon dow   command
#0 2,5,8,11,14,20 * * * /home/<user>/restart-terminated-instances.sh &>> /var/log/restart_gce_align_nodes/cron.log

instance_list=$(gcloud compute instances list --regexp "^align-([0-9]*)" | grep TERMINATED)

echo "[$(date)] Restarting "$(wc -l <<< "$instance_list")" nodes ..."

grep us-central1-b <<< "$instance_list" | awk '{print $1}' | xargs --no-run-if-empty gcloud compute instances start --zone us-central1-b
grep us-central1-c <<< "$instance_list" | awk '{print $1}' | xargs --no-run-if-empty gcloud compute instances start --zone us-central1-c
grep us-central1-f <<< "$instance_list" | awk '{print $1}' | xargs --no-run-if-empty gcloud compute instances start --zone us-central1-f

echo "[$(date)] Done."