#!/bin/bash
#SBATCH --job-name=exatractSNVsGnomad
#SBATCH --partition=cnu
#SBATCH --mem=100G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/data/encode/public/cosmicGnomad_20230215
#SBATCH --account=sscm013903

# Load modules for bed tools intersects
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc

get_cosmic_gnomad_overlaps() {
    # Get cosmic variants only, but keep original file to query other information for variants in the future
    # Create a 3000bp window around the variant to query against gnomAD data
    # Only keeps first instance if variants are duplicated
    cat $1 | awk -F"\t" '{print $1 "\t" $2-1000 "\t" $3+1000 "\t" $4 "\t" $5 "\t" $2 "\t" $NF}' | sed 1d | awk '$1 != "X" {print $0}' | awk '$1 != "Y" {print $0}' | awk '$1 != "MT" {print $0}' | awk '!visited[$0]++' > cosmic_snvs.bed
    sed -i -e 's/^/chr/' cosmic_snvs.bed

    # Sort files for intersecting
    bedtools sort -i cosmic_snvs.bed > cosmic.tmp
    mv cosmic.tmp cosmic.sorted.bed 
    bedtools sort -i $2 > gnomad.sorted.bed

    bedtools intersect -wa -wb -a cosmic.sorted.bed -b gnomad.sorted.bed > overlaps.bed

    # separating cosmic/gnomad variants, and removing duplicated gnomad variants from intersect command
    cat overlaps.bed | awk '{print $8 "\t" $9 "\t" $9 "\t" $11 "\t" $12 "\t" $7}' | awk '! a[$1, $2, $3, $4, $5]++' | sort > gnomad.bed
    cat overlaps.bed | awk '{print $1 "\t" $6 "\t" $6 "\t" $4 "\t" $5 "\t" $7}' | awk '! a[$0]++' | sort > cosmic.bed

    # convert to conservation format in the new directory
    cat cosmic.bed | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' | sort > cosmic.sorted.bed
    cat gnomad.bed | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' | sort > gnomad.sorted.bed

    # Remove variants that are common between both files
    # Add 0 and 1 for driver status
    comm -23 cosmic.sorted.bed gnomad.sorted.bed | awk '{print $0 "\t" 1}' > cosmic.tmp
    comm -23 gnomad.sorted.bed cosmic.sorted.bed | awk '{print $0 "\t" 0}' > gnomad.tmp
    mv cosmic.tmp cosmic1.sorted.bed
    mv gnomad.tmp gnomad1.sorted.bed

    # Merge the variant files
    cat cosmic1.sorted.bed gnomad1.sorted.bed | bedtools sort -i > /bp1/mrcieu1/data/encode/public/CanDrivR/training/$3/cosmicGnomadVariants_$3.bed

    ################ Reformat variants for conservation features ###################

    # convert to conservation format in the new directory
    cat cosmic.bed | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' | sort > cosmic.sorted.bed
    cat gnomad.bed | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' | sort > gnomad.sorted.bed

    # Remove variants that are common between both files
    # Add 0 and 1 for driver status
    comm -23 cosmic.sorted.bed gnomad.sorted.bed | awk '{print $0 "\t" 1}' > cosmic.tmp
    comm -23 gnomad.sorted.bed cosmic.sorted.bed | awk '{print $0 "\t" 0}' > gnomad.tmp
    mv cosmic.tmp cosmic1.sorted.bed
    mv gnomad.tmp gnomad1.sorted.bed

    # Merge the variant files
    cat cosmic1.sorted.bed gnomad1.sorted.bed | bedtools sort -i > /bp1/mrcieu1/data/encode/public/CanDrivR/training/$3/cosmicGnomadVariantsReformatted_$3.bed
    rm cosmic1.sorted.bed gnomad1.sorted.bed overlaps.bed cosmic_snvs.bed

}

#get_cosmic_gnomad_overlaps /bp1/mrcieu1/data/encode/public/CanDrivR/training/coding/cosmic/CosmicMutantExport.bed /bp1/mrcieu1/data/encode/public/CanDrivR/training/coding/gnomad/gnomad_exomes_snvs_0.05.sorted.bed coding
get_cosmic_gnomad_overlaps /bp1/mrcieu1/data/encode/public/CanDrivR/training/non-coding/cosmic/CosmicNCV.bed /bp1/mrcieu1/data/encode/public/CanDrivR/training/non-coding/gnomad/gnomad_genomes_snvs_0.05.sorted.bed non-coding
