#!/bin/bash
#SBATCH --job-name=getFeatureData
#SBATCH --partition=gpu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903
#SBATCH --array=1-35

# Load bedtools
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc

# Load the directory and file that you wish to query
featureDir="/bp1/mrcieu1/data/ucsc/public/ConsAll/"

# Get all of the available conservation feature datasets for querying
data=($featureDir*)

# Submit a different array job for each feature dataset
feature=$(echo ${data[@]} | cut --delimiter " " --fields ${SLURM_ARRAY_TASK_ID})

# Get the feature without the directory path
arrIN=(${feature//// })

# Specify output directory
variantDir="/bp1/mrcieu1/data/encode/public/1000G_ICGC/"
outputDir="/bp1/mrcieu1/data/encode/public/1000G_ICGC/features/"

# Load the positive and negative datasets
#posData=${variantDir}MSK_TCGA_ICGC_Variants.sorted.bed
#negData=${variantDir}1000GVariants.sorted.bed

variants=${variantDir}MSK_TCGA_ICGC_1000G.bed

# find intersects between cosmic/gnomad data and conservation scores
#bedtools intersect -wa -wb -a $feature -b $negData -sorted > $variantDir/neg.out.bed
#bedtools intersect -wa -wb -a $feature -b $posData -sorted > $variantDir/pos.out.bed

bedtools intersect -wa -wb -a $feature -b $variants -sorted > $outputDir/cons.out.bed

# reformat output of conservation intersect
#cat $outputDir/neg.out.bed | awk '{print $5 "\t" $6 "\t" $8 "\t" $9 "\t" $4}' > $outputDir/neg_${arrIN[6]}
#cat $outputDir/pos.out.bed | awk '{print $5 "\t" $6 "\t" $8 "\t" $9 "\t" $4}' > $outputDir/pos_${arrIN[6]}

cat $outputDir/cons.out.bed | awk '{print $5 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $4}' > $outputDir/${arrIN[6]}
