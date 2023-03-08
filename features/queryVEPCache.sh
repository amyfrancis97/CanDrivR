#!/bin/bash
#SBATCH --job-name=queryVEPCache
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

# Split the variant file into multiple files of ~100,000 lines each
# Must be broken down as query is too large for VEP
#mkdir /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit
#cd /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit
#cp /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/cosmicGnomadVariants.bed /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit
#split -l 100000 cosmicGnomadVariants.bed  cosmicGnomadVariants

# To query the cache, you must be in the emsembl-vep directory
# You must also have samtools/tabix/bedtools loaded in the environment
#cd ~/ensembl-vep

# Query VEP cache for consequence features
module load apps/bedops apps/bedtools
module load apps/bayestraits apps/bayestraits apps/bcftools apps/samtools apps/tabix lib/htslib

#files=(/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/*)

#for file in ${files[@]}; do

# Query VEP cache for consequence features
#./vep -i $file --cache --force_overwrite --vcf --fields "Consequence" -o ${file}_variant_effect_output_conseq.txt

# Query VEP cache for AA features
#./vep -i $file --cache --force_overwrite --vcf --fields "Amino_acids" -o ${file}_variant_effect_output_AA.txt

# Query VEP cache for distance features
#./vep -i $file --cache --force_overwrite --vcf --fields "DISTANCE" -o ${file}_variant_effect_output_distance.txt;
#done

#cat /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/*variant_effect_output_conseq-1.txt > /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/variant_effect_output_conseq.txt
#cat /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/*variant_effect_output_AA.txt > /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/variant_effect_output_AA.txt
#cat /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/*variant_effect_output_distance.txt > /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/variant_effect_output_distance.txt

rm /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/vepFileSplit/*cosmicGnomadVariant*




# Navigate back to script directory
#cd /bp1/mrcieu1/users/uw20204/paper1/features

# Reformat vep features to one-hot-encoding
#python reformatVEP.py
