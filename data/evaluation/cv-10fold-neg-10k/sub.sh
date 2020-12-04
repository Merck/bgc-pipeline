#!/bin/bash

split=$1

if [ -z "$split" ]; then
  echo "Usage: $0 split"
  exit 1
fi

if [ -f ${split}.job ]; then
   qstat -w $(cat ${split}.job) 2>&1 && echo "Job exists in queue, skipping..." && exit 1
fi

if [ -f ${split}.pickle ]; then
   echo "${split}.pickle exists, delete it to continue."
   exit 1
fi

rm ${split}.out 2>/dev/null
rm ${split}.err 2>/dev/null
rm ${split}.job 2>/dev/null
rm ${split}.test.csv 2>/dev/null

echo Submitting $(pwd)/../train_split.sh ${split}

wd=$(pwd)
jobname=${split}_$(basename ${wd})
qsub -V -q huge -l select=1:ncpus=8 -N ${jobname} -o ${split}.out -e ${split}.err -- $(pwd)/../train_split.sh $split | tee ${split}.job

exit $?
