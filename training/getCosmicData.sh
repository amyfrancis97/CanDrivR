#!/bin/bash
#SBATCH --job-name=getGnomad
#SBATCH --partition=cnu
#SBATCH --time=14-00:00:0
#SBATCH --mem=100G
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/CanDrivR
#SBATCH --account=sscm013903

module load apps/bedops apps/bedtools

############### Get Cosmic Coding Data ############
# Get the cosmic dataset
#wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicMutantExport.tsv.gz

# unzip
#gunzip CosmicMutantExport.tsv.gz

############### Get Cosmic Non-Coding Data ##########
#wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicNCV.tsv.gz

# unzip
#gunzip CosmicNCV.tsv.gz

filter_cosmic() {

echo $1

# Convert from TSV to CSV
tr '\t' ',' < $1.tsv > $1.csv 

# Remove white spaces in column names
# Select only SNVs (those with arrow symbol)
sed '1s/ /_/g'  $1.csv | awk -F, -v var=$2 '$var ~ '/'>'/'' > $1.tmp 
mv $1.tmp $1.csv

# Gets recurrence of each SNV (how many different sample ID contain the variant)
python /bp1/mrcieu1/users/uw20204/paper1/training/getReccurence.py $1 $2

}
# Execute function for cosmic coding and non-coding
# Number represents HGVSG column in dataframe

cd /bp1/mrcieu1/data/encode/public/CanDrivR/training/non-coding/cosmic
filter_cosmic CosmicNCV 24

cd /bp1/mrcieu1/data/encode/public/CanDrivR/training/coding/cosmic
filter_cosmic CosmicMutantExport 37
