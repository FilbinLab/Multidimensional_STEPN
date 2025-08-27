#!/bin/bash

if [[ ${#} -ne 1 ]]
  then
  echo "usage: source 3_NicheAnalysis.sh SampleName"
  fi
  
  SampleName=${1}
  
  script_dir=/n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions
  log_dir=/n/scratch/users/s/sad167/EPN/Xenium/logs/3_NicheAnalysis
  
  FILE=${log_dir}/NicheAnalysis-${SampleName}.out
  if [ ! -f "$FILE" ]; then
  sbatch -o ${log_dir}/NicheAnalysis-${SampleName}.out -e ${log_dir}/NicheAnalysiss-${SampleName}.err -J NicheAnalysis-${SampleName} ${script_dir}/3b_NicheAnalysis.sbatch ${script_dir} ${SampleName}
  fi
  