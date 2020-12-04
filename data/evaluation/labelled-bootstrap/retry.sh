#!/bin/bash

models=$@

echo Retrying models: ${models}

submitted=0
for file in ${models}; do 
  for index in 0 1 2 3 4; do 
    echo "-----"
    echo $file $index
    cd $file; 
    ../sub.sh $index && submitted=$((submitted + 1)); 
    cd ..; 
  done; 
done

echo ""
echo "Submitted ${submitted} jobs."
