#!/bin/bash
#SBATCH --job-name=getCosmicData
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

# Filter to only include substitutions
cat CosmicMutantExport.tsv| awk -F '\t' '$22 ~ /Substitution/' > CosmicMutantExportSubst.tsv

# Convert from TSV to CSV
tr '\t' ',' < CosmicMutantExport.tsv > CosmicMutantExport.csv 

# Remove white spaces in column names
sed '1s/ /_/g'  CosmicMutantExport.csv > CosmicMutantExportNoSpaces.csv 
mv CosmicMutantExportNoSpaces.csv CosmicMutantExport.csv

# Filter to only include substitutions
awk -F, 'NR==1 || $22 ~ /Substitution/' CosmicMutantExport.csv > CosmicMutantExportSubs.csv

# Partly extract the chromosome, position, ref and alternative alleles
cut -d',' -f37 CosmicMutantExportSubs.csv |awk -F  '[:,., >, -]' 'BEGIN{printf "CHROM\tX\tPOSREF\tALT\n"} {printf ("%s\t %s\t %s\t %s\n" ,$1, $2, $3, $4)}' > CosmicMutantExportSubs2.csv

# Further splits the pos and ref alleles by string position
# Then counts the reccurence of the variants
# Exports new file as BED "CosmicMutantExport.bed"
python getReccurence.py

# From the cosmic file, prints out: chromosome, position, position, reference allele, alternate allele
cat /user/home/uw20204/data/encode/public/cosmic/CosmicMutantExport.bed | awk '{print $1 "\t" $2 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > /user/home/uw20204/data/encode/public/cosmic/cosmic.bed
