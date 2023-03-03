#!/bin/bash
#SBATCH --job-name=getGnomad
#SBATCH --partition=cnu
#SBATCH --time=14-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215
#SBATCH --account=sscm013903

module load apps/bedops apps/bedtools

###############################################################	Cosmic Data #############
# Get the cosmic dataset
#wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicMutantExport.tsv.gz

# unzip
#gunzip CosmicMutantExport.tsv.gz

# Convert from TSV to CSV
#tr '\t' ',' < CosmicMutantExport.tsv > CosmicMutantExport.csv 

# Remove white spaces in column names
# Select only SNVs (those with arrow symbol)
#sed '1s/ /_/g'  CosmicMutantExport.csv | awk -F, '$37 ~ '/'>'/'' > CosmicMutantExport.tmp 
#mv CosmicMutantExportNoSpaces.csv CosmicMutantExport.csv

# Gets recurrence of each SNV (how many different sample ID contain the variant)
python /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/finalV2Scripts/getReccurence.py

