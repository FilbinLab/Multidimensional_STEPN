#!/bin/bash

if [[ ${#} -ne 1 ]]
  then
  echo "usage: source 2_assign_programs.sh SampleName"
  fi
  
  SampleName=${1}
  
  script_dir=/n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions
  log_dir=/n/scratch/users/s/sad167/EPN/Xenium/logs/2_assign_programs
  
  FILE=${log_dir}/AssignPrograms-${SampleName}.out
  if [ ! -f "$FILE" ]; then
  sbatch -o ${log_dir}/AssignPrograms-${SampleName}.out -e ${log_dir}/AssignPrograms-${SampleName}.err -J AssignPrograms-${SampleName} ${script_dir}/2b_assign_programs_automatic.sbatch ${script_dir} ${SampleName}
  fi
  