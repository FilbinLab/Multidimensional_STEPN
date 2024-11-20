#!/bin/bash

if [[ ${#} -ne 1 ]]
  then
  echo "usage: source 1_preprocessing.sh SampleName"
  fi
  
  SampleName=${1}
  
  script_dir=/n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions
  log_dir=/n/scratch/users/s/sad167/EPN/Xenium/logs/1_preprocessing
  
  FILE=${log_dir}/Prep-${SampleName}.out
  if [ ! -f "$FILE" ]; then
  sbatch -o ${log_dir}/Prep-${SampleName}.out -e ${log_dir}/Prep-${SampleName}.err -J Prep-${SampleName} ${script_dir}/1b_preprocessing_v2.sbatch ${script_dir} ${SampleName}
  fi
  