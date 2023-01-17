#!/bin/bash
#SBATCH --job-name=getCosmicGnomadOverlaps
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/user/home/uw20204/data/encode/public/gnomad

# Load bedtools
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc


# Sorts both the cosmic and gnomad datasets
sort /user/home/uw20204/data/encode/public/cosmic/cosmic.bed > /user/home/uw20204/data/encode/public/cosmic/cosmic.sorted.bed
sort cosmicGnomadControlsAndNonCancer/gnomad.bed > cosmicGnomadControlsAndNonCancer/gnomad.sorted.bed

# Test to see if there is any overlap between the gnomad & cosmic variants
# If there are, then remove the variants from cosmic as they are likely benign
comm -23 /user/home/uw20204/data/encode/public/cosmic/cosmic.sorted.bed cosmicGnomadControlsAndNonCancer/gnomad.sorted.bed > cosmicGnomadControlsAndNonCancer/cosmicGnomadExcl.bed

# Reformat cosmic variant file by making the start position 1000bp lower, and the end position 1000bp higher
# This creates a new dataset for querying againt the gnomad dataset
# Ensures gnomad variants overlap within 1000 bp of a cosmic variant
cat cosmicGnomadControlsAndNonCancer/cosmicGnomadExcl.bed | awk '{print $1 "\t" $2 - 1000 "\t" $2 + 1000 "\t" $4 "\t" $5 "\t" $2 "\t" $6}' > cosmicGnomadControlsAndNonCancer/cosmic.bed

# Intersect the new Cosmic dataset with the gnomad dataset, to pull out gnomad variants that overlap within 1000bp
bedtools intersect -wa -wb -a cosmicGnomadControlsAndNonCancer/cosmic.bed -b cosmicGnomadControlsAndNonCancer/gnomad.sorted.bed > cosmicGnomadControlsAndNonCancer/overlaps.bed

# Separate cosmic/gnomad variants from the overlaps file, and remove duplicated gnomad variants from intersect command 
#(since a single gnomad variant may fall within several cosmic variants 1000bp regions if the cosmic variants are close together)
cat cosmicGnomadControlsAndNonCancer/overlaps.bed | awk '{print $8 "\t" $9 "\t" $9 "\t" $11 "\t" $12 "\t" $7}' | awk '! a[$1, $2, $3, $4, $5]++' | sort > cosmicGnomadControlsAndNonCancer/gnomad.bed
cat cosmicGnomadControlsAndNonCancer/overlaps.bed | awk '{print $1 "\t" $6 "\t" $6 "\t" $4 "\t" $5 "\t" $7}' | awk '! a[$0]++' | sort > cosmicGnomadControlsAndNonCancer/cosmic.bed

# Convert to conservation format in the new directory (in the form pos, pos + 1)
mkdir cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/
cat cosmicGnomadControlsAndNonCancer/cosmic.bed | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $4 "\t" $5 "\t" $6}' > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmicQuery.bed
cat cosmicGnomadControlsAndNonCancer/gnomad.bed | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $4 "\t" $5 "\t" $6}' > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomadQuery.bed

# Sort formatted variants
bedtools sort -i cosmicGnomadControlsAndNonCancer/cosmic.bed > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomad.sorted.bed
bedtools sort -i  cosmicGnomadControlsAndNonCancer/gnomad.bed > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmic.sorted.bed
bedtools sort -i  cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomadQuery.bed > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomadQuery.sorted.bed
bedtools sort -i  cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmicQuery.bed > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmicQuery.sorted.bed

# Final datasets are: gnomad.sorted.bed, cosmic.sorted.bed, gnomadQuery.sorted.bed, cosmicQuery.sorted.bed
