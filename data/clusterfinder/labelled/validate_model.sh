#!/bin/bash

set -e

model=$1
color=$2
suffix=$3

filename=$(basename "$model")

input=../CF_labelled_contigs_domains${suffix}.csv
csv_output=${filename/.pickle/}${suffix}.csv
png_output=${filename/.pickle/}${suffix}.png

echo Validating ${model} with color ${color}, suffix "${suffix}"
echo ${input}
echo ${csv_output}
echo ${png_output}

dvc run -d ${model} -d ${input} -o ${csv_output} \
  python ../../../../bgc_detection/run_prediction.py -m ${model} -o ${csv_output} -e 0.01 ${input}

dvc run -d ${csv_output} -O ${png_output} \
  "python ../../../../bgc_detection/evaluation/prediction_roc.py --prediction ${csv_output} \
  --name ${filename/.pickle/} --color ${color} --title 'Labelled contig per-domain ROC' -o ${png_output}"