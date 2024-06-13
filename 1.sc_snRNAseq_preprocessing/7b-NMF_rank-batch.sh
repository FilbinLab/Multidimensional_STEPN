#!/bin/bash

if [[ ${#} -ne 1 ]]
then
        echo "usage: source NMF-batch.sh sample_name"
fi

sample_name=${1}

script_dir=/n/scratch/users/s/sad167/EPN/scRNAseq/scripts
log_dir=/n/scratch/users/s/sad167/EPN/scRNAseq/logs/NMF

FILE=${log_dir}/rank6_${sample_name}.out
if [ ! -f "$FILE" ]; then
    sbatch -o ${log_dir}/rank6_${sample_name}.out -e ${log_dir}/rank6_${sample_name}.err -J NMF_rank6-${sample_name} ${script_dir}/7a-NMF_rank.sbatch ${script_dir} ${sample_name}
fi
