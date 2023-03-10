#!/bin/bash
#SBATCH --job-name=getGC_CpG
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

python getGC_GpG.py $1
