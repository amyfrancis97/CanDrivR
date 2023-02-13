#!/bin/bash
#SBATCH --job-name=getAAProperties
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

## Download aaSEA AAindex AA property data
wget https://github.com/cran/aaSEA/raw/master/data/AAindex.rda

Rscript AASubMatrices.R

Rscript extractAAproperties.R
