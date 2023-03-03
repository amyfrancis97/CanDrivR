#!/bin/bash
#SBATCH --job-name=exatractSNVsGnomad
#SBATCH --partition=cnu
#SBATCH --mem=90G
#SBATCH --time=4-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215
#SBATCH --account=sscm013903

module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc

################ EXOME #################

# Get the Gnomad exome dataset (GRCh38 liftover)
#wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
#mv gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz
#gunzip gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz

################### GENOME ################

# Get the Gnomad genome dataset (GRCh38 liftover)
#wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
#mv gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz
#gunzip gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz

##############################################################################################
# The following function pulls out variants present in controls or non-cancer data
# Then filters data to only include SNVs
# Then generates a file that only includes chrom, start, end, ref_allele, alt_allele

get_variants() {

arrIN=(${1//./ })

# Remove first line/ header of the gnomad file
sed '/^#/d' $1 > ${arrIN[1]}_${2}_gnomad.bed

# pulls out all controls as well as non-cancer samples
# Filters to SNVs
cat ${arrIN[1]}_${2}_gnomad.bed | awk '($8 ~ '/non_cancer/') || ($8 ~ '/controls/') {print $0}' | awk 'length($4) == 1 && length($5) == 1 { print }' > ${arrIN[1]}_${2}_gnomadAlleleFreqIncl.txt

# get variants where allele frequency is >0.05 in non-cancer cohort
python /bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215/finalV2Scripts/getAlleleFrequenciesFunc.py $2 ${arrIN[1]}_

cat ${arrIN[1]}_${2}_gnomad_AF_filtered.out | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}'  | awk 'length($4)== 1 && length($5) == 1'| bedtools sort -i  > gnomad_${arrIN[1]}_snvs_${2}.sorted.bed
rm ${arrIN[1]}_${2}_gnomad_AF_filtered.out ${arrIN[1]}_${2}_gnomadAlleleFreqIncl.txt ${arrIN[1]}_${2}_gnomad.bed ${arrIN[1]}_${2}_gnomad.bed

}

executeFunc() {

# Execute above function for both exome and genome data in gnomAD
#get_variants gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf $1
get_variants gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf $1

# Concatenate the gnomad data files and only keep first occurence if variants are duplicated in both files
# Remove "chr" string from chrom column
cat gnomad_exomes_snvs_${1}.sorted.bed gnomad_genomes_snvs_${1}.sorted.bed | awk '!visited[$0]++' | awk '{ gsub(/chr/,"", $1); print } ' > gnomad_snvs_$1.bed
}

# Execute function for different specified gnomAD allele frequencies
# Final gnomAD file will be gnomad_snvs_${AF}.bed
AF=0.05
executeFunc $AF
