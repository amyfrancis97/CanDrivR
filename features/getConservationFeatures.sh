#!/bin/bash
#SBATCH --job-name=getFeatureData
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=5-00:00:0
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
variantDir="/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210/"
outputDir="/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210/features/"

# Load the positive and negative datasets
variants=${variantDir}cosmicGnomadVariantsReformatted.bed

# find intersects between cosmic/gnomad data and conservation scores
bedtools intersect -wa -wb -a $feature -b ${variantDir}variants_reformatted.bed -sorted > $outputDir/${arrIN[6]}_cons.out.bed

# reformat output of conservation intersect
cat $outputDir/${arrIN[6]}_cons.out.bed | awk '{print $5"\t"$6"\t"$8"\t"$9"\t"$10"\t"$11"\t"$4}' > $outputDir/${arrIN[6]}

rm $outputDir/${arrIN[6]}_cons.out.bed
