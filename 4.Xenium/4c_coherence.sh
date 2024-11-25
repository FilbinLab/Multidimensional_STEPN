#!/bin/bash

if [[ ${#} -ne 1 ]]
  then
  echo "usage: source 4_coherence.sh SampleName"
  fi
  
  SampleName=${1}
  
  script_dir=/n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions
  log_dir=/n/scratch/users/s/sad167/EPN/Xenium/logs/4_coherence
  
  FILE=${log_dir}/Coherence-${SampleName}.out
  if [ ! -f "$FILE" ]; then
  sbatch -o ${log_dir}/Coherence-${SampleName}.out -e ${log_dir}/Coherence-${SampleName}.err -J Coherence-${SampleName} ${script_dir}/4b_coherence.sbatch ${script_dir} ${SampleName}
  fi
  