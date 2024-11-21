#!/bin/bash
#SBATCH -c 6                              # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
#SBATCH -t 0-01:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=180G                          # Memory total in MB (for all cores)
#SBATCH --job-name %j 		# Job name
#SBATCH -o %j.out		# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written

# Load required modules
module load gcc/9.2.0 cmake/3.22.2 R/4.4.0 python/3.10.11 gdal/3.1.4 udunits/2.2.28 geos/3.10.2 imageMagick/7.1.0 texlive/2007 git/2.35.3 hdf5/1.14.0

# Run the R script
Rscript ${1}/2_assign_programs_automatic.R ${2}
