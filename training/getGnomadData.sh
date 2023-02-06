#!/bin/bash
#SBATCH --job-name=getGnomadData
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/user/home/uw20204/data/encode/public/gnomad

# Load bedtools
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc

# Get the Gnomad exome dataset (GRCh38 liftover)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
mv gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz 
gzip -d gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz 

# Remove first line/ header of the gnomad file
sed '900d' gnomad.exomes.r2.1.1.sites.vcf > gnomad.exomes.r2.1.1.sites.bed

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
cat cosmicGnomadControlsAndNonCancer/gnomadAlleleFreqIncl.txt | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}'  | awk 'length($4)== 1 && length($5) == 1' > cosmicGnomadControlsAndNonCancer/gnomad.bed
