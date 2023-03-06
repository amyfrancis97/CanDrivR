#!/bin/bash
#SBATCH --job-name=queryVEPCache
#SBATCH --partition=cnu
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

./vep -i /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/cosmicGnomadVariants.bed --cache --force_overwrite --vcf --fields "Consequence" -o "variant_effect_output_conseq.txt"

# Query VEP cache for AA features
./vep -i /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/cosmicGnomadVariants.bed --cache --force_overwrite --vcf --fields "Amino_acids" -o "variant_effect_output_AA.txt"

# Query VEP cache for distance features
./vep -i /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/cosmicGnomadVariants.bed --cache --force_overwrite --vcf --fields "DISTANCE" -o "variant_effect_output_distance.txt"

# Navigate back to script directory
cd /bp1/mrcieu1/users/uw20204/paper1/features

# Reformat vep features to one-hot-encoding
python reformatVEP.py
