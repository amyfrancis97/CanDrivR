#!/bin/bash
#SBATCH --job-name=getDNAShape
#SBATCH --partition=cnu
#SBATCH --mem=100G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

module load lang/r
Rscript DNAShape.R $1
