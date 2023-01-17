#!/bin/bash
#SBATCH --job-name=exatractSNVsGnomad
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/user/home/uw20204/data/encode/public/gnomad

# Load bedtools
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc

# Get the cosmic dataset
wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicMutantExport.tsv.gz

# unzip
gunzip CosmicMutantExport.tsv.gz

# Convert from TSV to CSV
tr '\t' ',' < CosmicMutantExport.tsv > CosmicMutantExport.csv 

# Remove white spaces in column names
sed '1s/ /_/g'  CosmicMutantExport.csv > CosmicMutantExportNoSpaces.csv 
mv CosmicMutantExportNoSpaces.csv CosmicMutantExport.csv

# Filter to only include substitutions
awk -F, '$22 ~ /Substitution/' CosmicMutantExport.csv > CosmicMutantExportSubs.csv

# Get the Gnomad exome dataset (GRCh38 liftover)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz

# Remove first line/ header of the gnomad file
sed '1d' gnomad.exomes.r2.1.1.sites.liftover_grch38.bed > gnomad.exomes.r2.1.1.sites.liftover_grch38.2.bed

# Pulls out all rows where there is a control as well as non-cancer column
cat gnomad.exomes.r2.1.1.sites.liftover_grch38.2.bed | awk '($8 ~ '/non_cancer_AF=/') || ($8 ~ '/controls_AF=/') {print $0}' > gnomadAlleleFreqIncl.txt

# Get variants where allele frequency is >0.05 in non-cancer/control cohort cohort
# Filters to more common variants to have a greater chance that these will indeed be neutral
python getAlleleFrequencies.py

# test to see if there is any overlap between the gnomad & cosmic variants
# if there are, then remove from cosmic as they are likely benign

# From the new gnomad variant file that filters to allele frequency > 5%
# Filters to single nucleotide variants only
# Pulls out chromosome, position, position, reference allele, alternate allele
cat cosmicGnomadControlsAndNonCancer/gnomadAlleleFreqIncl.txt | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}'  | awk 'length($4)== 1 && length($5) == 1' > cosmicGnomadControlsAndNonCancer/gnomadMatch.bed

# From the cosmic file, prints out: chromosome, position, position, reference allele, alternate allele
cat /user/home/uw20204/data/encode/public/cosmic/CosmicMutantExportFilteredSubst.bed | awk '{print $1 "\t" $2 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > /user/home/uw20204/data/encode/public/cosmic/cosmicMatch.bed

# Sorts both the cosmic and gnomad datasets
sort /user/home/uw20204/data/encode/public/cosmic/cosmicMatch.bed > /user/home/uw20204/data/encode/public/cosmic/cosmicMatch.sorted.bed
sort cosmicGnomadControlsAndNonCancer/gnomadMatch.bed > cosmicGnomadControlsAndNonCancer/gnomadMatch.sorted.bed

# Test to see if there is any overlap between the gnomad & cosmic variants
# If there are, then remove the variants from cosmic as they are likely benign
comm -23 /user/home/uw20204/data/encode/public/cosmic/cosmicMatch.sorted.bed cosmicGnomadControlsAndNonCancer/gnomadMatch.sorted.bed > cosmicGnomadControlsAndNonCancer/cosmicGnomadExcl.bed

# Reformat cosmic variant file by making the start position 1000bp lower, and the end position 1000bp higher
# This creates a new dataset for querying againt the gnomad dataset
# Ensures gnomad variants overlap within 1000 bp of a cosmic variant
cat cosmicGnomadControlsAndNonCancer/cosmicGnomadExcl.bed | awk '{print $1 "\t" $2 - 1000 "\t" $2 + 1000 "\t" $4 "\t" $5 "\t" $2 "\t" $6}' > cosmicGnomadControlsAndNonCancer/cosmicQuery.bed

# Intersect the new Cosmic dataset with the gnomad dataset, to pull out gnomad variants that overlap within 1000bp
bedtools intersect -wa -wb -a cosmicGnomadControlsAndNonCancer/cosmicQuery.bed -b cosmicGnomadControlsAndNonCancer/gnomadMatch.sorted.bed > cosmicGnomadControlsAndNonCancer/overlaps.bed

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

# Final datasets are: gnomadQuery.bed and cosmicQuery.bed 
