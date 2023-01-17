#!/bin/bash
#SBATCH --job-name=exatractSNVsGnomad
#SBATCH --partition=cnu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/user/home/uw20204/data/encode/public/gnomad

# Load bedtools
module load apps/bedops/2.4.38 apps/bedtools/2.30.0
module load apps/bcftools apps/samtools/1.9 apps/tabix/0.2.5 lib/htslib/1.10.2-gcc
# The options for gnomad data can be:
# "cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomad.sorted.bed"
# "cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomadQuery.sorted.bed" - this is for conservation where pos, pos + 1
# The options for cosmic data can be:
# "cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmic.sorted.bed"
# "cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmicQuery.sorted.bed" - this is for conservation where pos, pos + 1

# Load the directory and file that you wish to query
dir="/user/home/uw20204/data/ucsc/public/k50.Umap.MultiTrackMappability/released/2021-11-01/"
file="k50.Umap.MultiTrackMappability.bedGraph"
cosmicData="cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmicQuery.sorted.bed"
gnomadData="cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomadQuery.sorted.bed"

# find intersects between cosmic/gnomad data and conservation scores
conservationDataNoSort=$($dir$file)

bedtools sort -i $conservationDataNoSort > sorted_${conservationDataNoSort}
#rm $conservationDataNoSort

conservationData=sorted_${conservationDataNoSort}

bedtools intersect -wa -wb -a $gnomadData -b $conservationData -sorted > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomad.out.bed
bedtools intersect -wa -wb -a $cosmicData -b $conservationData -sorted > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmic.out.bed

# reformat output of conservation intersect
cat cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomad.out.bed | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $10}' > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/gnomad$file
cat cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmic.out.bed | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $10}' > cosmicGnomadControlsAndNonCancer/queryPhyloPhastCons/cosmic$file
