#!/bin/bash
#SBATCH --job-name=getDinucleotides
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

# Get conservation features
#sbatch getConservationFeatures.sh

# Get VEP features
sbatch queryVEPCache.sh


Rscript dinucleotideProperties.R
