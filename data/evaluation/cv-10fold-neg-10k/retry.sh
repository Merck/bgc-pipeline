#!/bin/bash

models=$@

echo Retrying models: ${models}

submitted=0
for file in ${models}; do 
    for split in 1 2 3 4 5 6 7 8 9 10; do 
      echo "-----"
      echo ${file}/fold${split}
      cd $file; 
      ../sub.sh fold${split} && submitted=$((submitted + 1)); 
      cd ..; 
    done
done

echo ""
echo "Submitted ${submitted} jobs."
