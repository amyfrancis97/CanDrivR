#!/bin/bash
#SBATCH --job-name=getGnomad
#SBATCH --partition=cnu
#SBATCH --time=14-00:00:0
#SBATCH --mem=50G
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210
#SBATCH --account=sscm013903

module load apps/bedops apps/bedtools

################ EXOME #################

# Get the GnomAD exome dataset (GRCh38 liftover)
#wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
#mv gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz 
#gzip -d gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz 

################### GENOME ################

# Get the Gnomad genome dataset (GRCh38 liftover)
#wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
#mv gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz
#gzip -d gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz

# The following function pulls out variants present in controls or non-cancer data
# Then filters data to only include SNVs
# Then generates a file that only includes chrom, start, end, ref_allele, alt_allele

#get_variants() {
# Remove first line/ header of the gnomad file
#sed '900d' $1 > gnomad.bed

# Pulls out all rows where there is a control as well as non-cancer column
#cat gnomad.bed | awk '($8 ~ '/non_cancer_AF=/') || ($8 ~ '/controls_AF=/') {print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > gnomad2.bed

# Filter to only include SNVs and exlude chromosomes X & Y
#cat gnomad2.bed | awk '(length($4) == 1) && (length($5) == 1) { print }' | awk '$1 != "X" {print $0}' | awk '$1 != "Y" {print $0}' | awk '$1 != "chrX" {print $0}' | awk '$1 != "chrY" {print $0}'> gnomad3.bed

#cd /bp1/mrcieu1/users/uw20204/paper1/training/

# Filter gnomad variants with AF > 0.01 for non-cancer/control
#python getGnomadAF.py

#cd /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230210

#arrIN=(${1//./ })

# Get gnomad variants only, but keep original file to query other information for variants in the future
#cat gnomadGreater0.01AF.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6}'| sed 1d > gnomad_${arrIN[1]}_snvs.bed
#}

# Execute above function for both exome and genome data in gnomAD
#get_variants gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf
#get_variants gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf

# Concatenate the gnomad data files and only keep first occurence if variants are duplicated in both files
# Remove "chr" string from chrom column
cat gnomad_exomes_snvs.bed gnomad_genomes_snvs.bed | awk '!visited[$0]++' | awk '{ gsub(/chr/,"", $1); print } ' > gnomad_snvs.bed

bedtools sort -i gnomad_snvs.bed >gnomad_snvs.sorted.bed
