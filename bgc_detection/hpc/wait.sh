#!/bin/bash
# Wait for all jobs to finish based on passed Job IDs
# Job IDs are passed as arguments, they should be separated by newlines or spaces

# Fail on error
set -e

JOB_LIST_FILE=$1

if [ -z "JOB_LIST_FILE" ]; then
  echo "Usage: $0 JOB_LIST_FILE"
  exit 1
fi

# Get status of given job IDs
function jobList() {
    qstat $(cat ${JOB_LIST_FILE}) 2>/dev/null | head | grep "$USER"
}

if ! type qstat > /dev/null 2>&1; then
    echo "Command qstat not available, not waiting."
    exit 1
fi

echo "Waiting for all jobs to finish..."
sleep 1;
while jobList; do
    echo "Waiting for jobs to finish..."
    sleep 30;
done
echo "Done, sleeping for 1 minute to make sure all output is copied..."
sleep 60;
