#!/bin/bash
#SBATCH -c 8                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40G                          # Memory total in MB (for all cores)

module load gcc/6.2.0 samtools/1.3.1 hisat2/2.1.0 rsem/1.3.0 R/4.1.1

sh /n/scratch/users/s/sad167/EPN/scripts/preprocessing/main.sh $1 $2
