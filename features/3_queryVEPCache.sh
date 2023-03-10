#!/bin/bash
#SBATCH --job-name=queryVEPCache
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

dataset=$1
dir=/bp1/mrcieu1/data/encode/public/CanDrivR/training/${dataset}/
# Split the variant file into multiple files of ~100,000 lines each
# Must be broken down as query is too large for VEP
newDir=${dir}vepFileSplit.tmp/
mkdir ${newDir}
cd ${newDir}
cp ${dir}cosmicGnomadVariants_${dataset}.bed ${newDir}

echo ${newDir}
pwd

split -l 100000 cosmicGnomadVariants_${dataset}.bed  cosmicGnomadVariants_${dataset}

# To query the cache, you must be in the ensembl-vep directory
# You must also have samtools/tabix/bedtools loaded in the environment
cd ~/ensembl-vep

# Query VEP cache for consequence features
module load apps/bedops apps/bedtools
module load apps/bayestraits apps/bayestraits apps/bcftools apps/samtools apps/tabix lib/htslib

files=(${newDir}*)

for file in ${files[@]}; do

# Query VEP cache for consequence features
./vep -i $file --cache --force_overwrite --vcf --fields "Consequence" -o ${file}_variant_effect_output_conseq.txt

# Query VEP cache for AA features
./vep -i $file --cache --force_overwrite --vcf --fields "Amino_acids" -o ${file}_variant_effect_output_AA.txt

# Query VEP cache for distance features
./vep -i $file --cache --force_overwrite --vcf --fields "DISTANCE" -o ${file}_variant_effect_output_distance.txt;
done

# Concatenate all of the vep output files
cat ${newDir}*variant_effect_output_conseq.txt > ${dir}features/variant_effect_output_conseq.txt
cat ${newDir}*variant_effect_output_AA.txt > ${dir}features/variant_effect_output_AA.txt
cat ${newDir}*variant_effect_output_distance.txt > ${dir}features/variant_effect_output_distance.txt

# Remove the temporary directory
#rm -r ${newDir}

sed -i '/^#/d' ${dir}features/variant_effect_output_conseq.txt 
sed -i '/^#/d' ${dir}features/variant_effect_output_distance.txt
sed -i '/^#/d' ${dir}features/variant_effect_output_AA.txt

# Navigate back to script directory
cd /bp1/mrcieu1/users/uw20204/paper1/features

# Reformat vep features to one-hot-encoding
#python reformatVEP.py
