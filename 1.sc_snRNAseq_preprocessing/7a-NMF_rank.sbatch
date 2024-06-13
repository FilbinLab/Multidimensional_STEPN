#!/bin/bash

#SBATCH -c 8                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -p short		# partition name
#SBATCH -t 0-12:00 		# hours:minutes runlimit after which job will be killed
#SBATCH --mem 32G 		# amount of memory requested
#SBATCH --job-name %j 		# Job name
#SBATCH -o %j.out		# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written

# Load required modules
module load gcc/9.2.0 R/4.2.1 geos/3.10.2

# Run the R script
Rscript ${1}/7-NMF_rank-forO2.R ${2}
