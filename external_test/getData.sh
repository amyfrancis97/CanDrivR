for i in {1..22}; do
## Download datasets for 1000 G for specific chromosome
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

## Unzip file
gunzip 1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

## Use BCF tools to extract relevant info
bcftools query -f '%CHROM %POS %REF %ALT %AF\n' 1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf > chr${i}_1000G.bed

## Extract variants where AF > 5% & only single nucleotide substitutions
cat chr${i}_1000G.bed |awk '{ if ($5 > 0.05) { print } }' | awk '{if (length($3) == 1 && length($4) == 1) { print }}' > chr${i}_1000G_Great5PercAF.bed

## Remove unwanted files
rm 1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf chr${i}_1000G.bed
;
done

# Concatenate the files for all chromosomes
cat *_Great5PercAF.bed  > 1000G.bed
rm *_Great5PercAF.bed 

# Make a bed file for start (-1000bp) and end (+1000bp) to query new negative dataset for overlaps within region of 1000 bp
cat 1000G.bed | awk '{print $1"\t"$2-1000"\t"$2+1000"\t"$3"\t"$4"\t"$2}' | bedtools sort -i > 1000GplusMin1000.bed

##################################################
# Get TCGA and ICGC and MSK positive data
# download all datasets of interest using download manager
# Downloaded all ICGC datasets from cBioPortal
rm *.html

getPositiveVariants() {
    # Unzip
    for i in $(ls *$1*); do tar -zxvf $i ${i%???????}/data_mutations.txt; done 

    # find all mutation files and concatenate
    cat $(find /Users/uw20204/Downloads -name data_mutations.txt) > allMutations.txt

    rm -r *tcga*

    # Filter to only contain SNPs
    cat allMutations1.txt | awk '$11 ~ /SNP/ { print }' > allMutationsSNP.txt

    cat allMutationsSNP.txt | awk '{print $5"\t"$6"\t"$7"\t"$12"\t"$14}' > GRCh37.bed

    # Liftover variants from hg19 to GRCh38
    python liftover.py

    # Concatenate all files together, sort and remove duplicate variants.
    cat *GRCh38.bed* | bedtools sort -i | uniq -u > GRCh38_${1}.sorted.bed
}
getPositiveVariants tcga
getPositiveVariants icgc
getPositiveVariants msk


# Concatenate TCGA, ICGC, and MSK variants, removing duplicted variants
cat GRCh38_tcga.sorted.bed GRCh38_icgc.sorted.bed  GRCh38_msk.sorted.bed |bedtools sort -i | uniq -u > positiveVariants.sorted.bed

# Get intersects between positive and negative examples
bedtools intersect -wa -wb -a positiveVariants.sorted.bed -b 1000GplusMin1000.sorted.bed -sorted > overlaps.bed

cat overlaps.bed | awk '{print $6"\t"$11"\t"$11"\t"$9"\t"$10}' | sort | uniq -u | bedtools sort -i > 1000GVariants.bed
cat overlaps.bed | awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5}' | sort | uniq -u | bedtools sort -i > ICGC_TCGA_MSK_Variants.bed
