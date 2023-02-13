#!/bin/bash
#SBATCH --job-name=getGnomad
#SBATCH --partition=cnu
#SBATCH --time=14-00:00:0
#SBATCH --mem=50G
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210
#SBATCH --account=sscm013903

module load apps/bedops apps/bedtools

###############################################################	Cosmic Data #############
# Get the cosmic dataset
#wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicMutantExport.tsv.gz

# unzip
#gunzip CosmicMutantExport.tsv.gz

# Gets recurrence of each SNV (how many different sample ID contain the variant)
#Python getReccurence.py

#gunzip CosmicMutantExportNew.bed.gz

# Get cosmic variants only, but keep original file to query other information for variants in the future
# Create a 3000bp window around the variant to query against gnomAD data
# Only keeps first instance if variants are duplicated
cat CosmicMutantExportNew.bed | awk -F"\t" '{print $1 "\t" $2-3000 "\t" $3+3000 "\t" $4 "\t" $5 "\t" $2 "\t" $NF}' | sed 1d | awk '$1 != "X" {print $0}' | awk '$1 != "Y" {print $0}' | awk '$1 != "MT" {print $0}' | awk '!visited[$0]++' > cosmic_snvs.bed

