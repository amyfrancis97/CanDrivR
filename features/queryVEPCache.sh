#!/bin/bash
#SBATCH --job-name=queryVEPCache
#SBATCH --partition=gpu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

# To query the cache, you must be in the emsembl-vep directory
# You must also have samtools/tabix/bedtools loaded in the environment
cd ~/ensembl-vep

# Query VEP cache for consequence features
module load apps/bedops apps/bedtools
module load apps/bayestraits apps/bayestraits apps/bcftools apps/samtools apps/tabix lib/htslib

./vep -i /bp1/mrcieu1/data/encode/public/1000G_ICGC/MSK_TCGA_ICGC_1000G.bed --cache --force_overwrite --vcf --fields "Consequence"

# Query VEP cache for AA features
./vep -i /bp1/mrcieu1/data/encode/public/1000G_ICGC/MSK_TCGA_ICGC_1000G.bed --cache --force_overwrite --vcf --fields "Amino_acids" --output_file "variant_effect_output_AA.txt"

# Query VEP cache for distance features
./vep -i /bp1/mrcieu1/data/encode/public/1000G_ICGC/MSK_TCGA_ICGC_1000G.bed --cache --force_overwrite --vcf --fields "DISTANCE" --output_file "variant_effect_output_distance.txt"

# Navigate back to script directory
cd /bp1/mrcieu1/users/uw20204/paper1/features

# Reformat vep features to one-hot-encoding
python reformatVEP.py
