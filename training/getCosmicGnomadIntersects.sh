#!/bin/bash
#SBATCH --job-name=getGnomad
#SBATCH --partition=cnu
#SBATCH --time=14-00:00:0
#SBATCH --mem=50G
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210
#SBATCH --account=sscm013903

# Load bedtools
module load apps/bedops apps/bedtools

# Sorts both the cosmic and gnomad datasets
bedtools sort -i cosmic_snvs.bed > cosmic_snvs.sorted.bed
bedtools sort -i /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210/gnomad/gnomad_exomes_snvs.bed > /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210/gnomad/gnomad_exomes_snvs.sorted.bed


# Get intersects between gnomad data and cosmic data
bedtools intersect -wa -wb -a cosmic_snvs.sorted.bed -b gnomad/gnomad_exomes_snvs.sorted.bed -sorted > intersectsGnomadCosmic.bed

# Get cosmic variants
cat intersectsGnomadCosmic.bed | awk '{print $1 "\t" $6 "\t" $6 "\t" $4 "\t" $5 "\t" $7}' | sort | uniq | sort > cosmic_overlap_snvs.sorted.bed

# Get gnomad variants
cat intersectsGnomadCosmic.bed | awk '{print $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $7}' | sort | uniq | sort > gnomad_overlap_snvs.sorted.bed

# Test to see if there is any overlap between the gnomad & cosmic variants
# If there are, then remove the variants from both datasets
comm -1 -3 cosmic_overlap_snvs.sorted.bed gnomad_overlap_snvs.sorted.bed > gnomad_snvs.bed
comm -1 -3  gnomad_overlap_snvs.sorted.bed cosmic_overlap_snvs.sorted.bed > cosmic_snvs.bed

# If duplicated, select occurence where the R value is the highest
# This is because gnomad variants might overlap with multiple cosmic variants
# Want to find the highest recurrence
cat gnomad_snvs.bed | sort -k6,6rn | sort -uk1,1 -uk2,2 > gnomad_snvs2.bed
mv gnomad_snvs2.bed gnomad_snvs.bed

# Merge the variant files but label '0' for gnomad and '1' for cosmic
cat cosmic_snvs.bed | awk '{print $0 "\t" 1}' > cosmic.bed
cat gnomad_snvs.bed | awk '{print $0 "\t" 0}' > gnomad.bed
cat cosmic.bed gnomad.bed > cosmicGnomadVariants.bed
bedtools sort -i cosmicGnomadVariants.bed > cosmicGnomadVariants.sorted.bed

# Reformat variants for conservation match
sed -e 's/^/chr/' cosmicGnomadVariants.sorted.bed | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > cosmicGnomadVariantsReformatted.bed
