#!/bin/bash
#SBATCH -c 6                              # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
#SBATCH -t 1-00:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                           # Partition to run in
#SBATCH --mem=100G                          # Memory total in MB (for all cores)
#SBATCH -o /n/scratch/users/s/sad167/EPN/Xenium/logs/1_assign_programs.out
#SBATCH -e /n/scratch/users/s/sad167/EPN/Xenium/logs/1_assign_programs.err
#SBATCH -J Xenium-assign-programs-niches

module load gcc/9.2.0 cmake/3.22.2 R/4.2.1 gdal/3.1.4 udunits/2.2.28 geos/3.10.2

Rscript /n/scratch/users/s/sad167/EPN/Xenium/scripts/plots/1_assign_programs_Xenium.R