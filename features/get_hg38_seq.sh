#!/bin/bash
#SBATCH --job-name=getFeatureData
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=5-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

wget -O hg38_seq.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38_seq.fa.gz
