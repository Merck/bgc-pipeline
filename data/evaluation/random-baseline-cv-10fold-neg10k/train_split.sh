#!/bin/bash

split=$1
cuda_visible=$2

# Exit on error
set -e

cd $PBS_O_WORKDIR

echo USER $USER
echo WD $(pwd)

module unload cuda/8.0.44;
module load cuda/9.0;
module load pbs/13.1.2;
module add anaconda3;
echo "Activating conda BGC environment..."
source activate bgc;

export KERAS_BACKEND=tensorflow

if [ ! -z "$cuda_visible" ]; then
  nvidia-smi --query-gpu=name,index,utilization.gpu,utilization.memory --format=csv
fi

echo PYTHON $(which python) version $(python --version)

python ../../../../bgc_detection/run_training.py -c config.json -e 0.01 -o ${split}.pickle ../splits/${split}.train.csv

python ../../../../bgc_detection/run_prediction.py -m ${split}.pickle -o ${split}.test.csv -e 0.01 --whole ../splits/${split}.test.csv
