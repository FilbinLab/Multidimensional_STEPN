#!/bin/bash

#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
# -N 1 means all cores will be on the same node)
#SBATCH -p gpu
#SBATCH --gres=gpu:1		# partition name
#SBATCH -t 0-05:00 		# hours:minutes runlimit after which job will be killed
#SBATCH --mem 150G 		# amount of memory requested
#SBATCH --job-name Xenium_rename_idents 		# Job name
#SBATCH -o /n/scratch/users/s/sad167/EPN/Xenium/logs/sbatch_out.out		# File to which standard out will be written
#SBATCH -e /n/scratch/users/s/sad167/EPN/Xenium/logs/sbatch_errors.err 		# File to which standard err will be written

set -e

# Load required modules
module load gcc/9.2.0 R/4.2.1 udunits/2.2.28 geos/3.10.2 python/3.8.12

# Run the scripts
Rscript /n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/Xenium_rename_idents_20231020__200939__BT2126_BT1745.R
Rscript /n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/Xenium_rename_idents_20231102__215055__7EP1_7EP41_3EP8.R
Rscript /n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/Xenium_rename_idents_20231107__203958__BT1717_BT775.R
Rscript /n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/Xenium_rename_idents_20231109__203408__BT1804_BT2169.R
Rscript /n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/Xenium_rename_idents_20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41.R