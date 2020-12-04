#!/bin/bash
# Execute any command within the conda bgc environment

cd $PBS_O_WORKDIR

if [ -z "$PBS_O_WORKDIR" ]; then
  echo "Environment variable PBS_O_WORKDIR not present, are you executing from a HPC job?"
  exit 1
fi

module add anaconda3
source activate bgc

echo Dir: $(pwd)
echo Executing: $@

exec "$@"
