#!/bin/bash
# Wait for all jobs of current user to finish

# Fail on error
set -e

function jobList() {
    qstat 2>/dev/null | head | grep "$USER"
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
