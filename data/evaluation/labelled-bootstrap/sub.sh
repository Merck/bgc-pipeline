#!/bin/bash

index=$1

if [ -z "$index" ]; then
  echo "Usage: $0 index"
  exit 1
fi

if [ -f model${index}.job ]; then
   qstat -w $(cat model${index}.job) 2>&1 && echo "Job exists in queue, skipping..." && exit 1
fi

if [ -f model${index}.pickle ]; then
   echo "model${index}.pickle exists, delete it to continue."
   exit 1
fi

rm model${index}.out 2>/dev/null
rm model${index}.err 2>/dev/null
rm model${index}.job 2>/dev/null
rm -r model${index}.tensorboard 2>/dev/null

echo Submitting $(pwd)/train_split.sh ${index}

wd=$(pwd)
jobname=$(basename ${wd})_${index}
qsub -V -q huge -l select=1:ncpus=8 -N ${jobname} -o model${index}.out -e model${index}.err  -- $(pwd)/train_split.sh $index | tee model${index}.job

exit $?
