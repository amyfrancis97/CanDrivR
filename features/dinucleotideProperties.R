# Get DNA shapes of 10 BP regions overlapping with variant
.libPaths("/bp1/mrcieu1/users/uw20204/paper1/features/RpackageLib") 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

BiocManager::install("BSgenome", quietly= TRUE, force = TRUE, INSTALL_opts = '--no-lock')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38', quietly= TRUE, dependencies=TRUE, force = TRUE, INSTALL_opts = '--no-lock')
BiocManager::install('DNAshapeR', quietly= TRUE, force = TRUE, INSTALL_opts = '--no-lock')
install.packages('dplyr', quietly= TRUE,repos = "http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

library('BSgenome')
library('BSgenome.Hsapiens.UCSC.hg38')
library('DNAshapeR')
library('GenomicRanges')
library('dplyr')
library('stringr')
library(tidyverse)
source("config.R")

args <- commandArgs()

# Reads in variant file in the format: "chrom", "start", "end", "ref", "alt", "R", "driver_stat"
variant_type = args[6]
variants = paste(paste(paste(paste(dir, variant_type, sep = ""), "cosmicGnomadVariants_", sep = "/"), variant_type, sep = ""), ".bed", sep = "") 
featureOutputDir=paste(paste(dir, variant_type, sep = ""), "/features/", sep = "")

# Import variants for shaping
variants=read.table(variants, sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "R", "driver_stat")

# Get the desired base pair range for DNA shape
variants[2] = variants[2]-5
variants[3] = variants[3]+5

# Make a GRRanges object
variants = makeGRangesFromDataFrame(variants)

# Get the 10bp fasta for each variant
getFasta(variants, BSgenome = Hsapiens, width = 3, filename = "VariantDinucleotides.fa")

# Import variants for shaping
variants = paste(paste(paste(paste(dir, variant_type, sep = ""), "cosmicGnomadVariants_", sep = "/"), variant_type, sep = ""), ".bed", sep = "")
variants=read.table(variants, sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "R", "driver_stat")
VariantDinucleotideWTSeq=read.table("VariantDinucleotides.fa")
toDelete <- seq(1, nrow(VariantDinucleotideWTSeq), 2)
variants = cbind(variants, VariantDinucleotideWTSeq[ -toDelete ,])

getMutantTrinucleotides = function(variantRow){
  mutantTrinucleotides = paste(substr(variants[variantRow, 8], 1, 1), variants[variantRow, 5], substr(variants[variantRow, 8], 3, 3), sep = "")
  return(mutantTrinucleotides)
}

# Carry out function to retrieve mutant trinucleotides for each variant
variantdfapply <- lapply(1:nrow(variants), getMutantTrinucleotides)

# Melt lists of variants into a dataframe
variantdf = do.call(rbind.data.frame, variantdfapply)
variants = cbind(variants, variantdf)
colnames(variants) = c("chrom", "start", "end", "ref_allele", "alt_allele", "R", "driver_stat", "WTtrinuc", "mutTrinuc")
print(head(variants))

# Read in dinucleotide properties
dinucleotideProperty=read.csv(dinucleotidePropertyTable)

# Get names of dinucleotide properties
dinucleotidePropertyNames = apply(dinucleotideProperty['PropertyName'],2,function(x)gsub('\\s+', '_',x))
print(length(dinucleotidePropertyNames))

# Function gets dinucleotideproperties for trinucleotide
getDinucleotidePropertyVector = function(variantRow){
  di1 = paste(substr(variants[variantRow, 'WTtrinuc'], 1, 1), substr(variants[variantRow, 'WTtrinuc'], 2, 2), sep = "") # position1-position2 of WT
  di2 = paste(substr(variants[variantRow, 'WTtrinuc'], 2, 2), substr(variants[variantRow, 'WTtrinuc'], 3, 3), sep = "") # position3-position4 of WT
  di3 = paste(substr(variants[variantRow, 'mutTrinuc'], 1, 1), substr(variants[variantRow, 'mutTrinuc'], 2, 2), sep = "") # position1-position2 of mutant
  di4 = paste(substr(variants[variantRow, 'mutTrinuc'], 2, 2), substr(variants[variantRow, 'mutTrinuc'], 3, 3), sep = "") # position3-position4 of mutant
  x = cbind(t(dinucleotideProperty[di1]), t(dinucleotideProperty[di2]))
  y = cbind(t(dinucleotideProperty[di3]), t(dinucleotideProperty[di4]))
  z = cbind(x, y)
  return(z)
}

# Carry out function to retrieve mutant dicleotide properties for each variant
# Four concatenated vectors
dinucleotides <- lapply(1:nrow(variants), getDinucleotidePropertyVector)

# Melt lists of variants into a dataframe
variantdf = do.call(rbind.data.frame, dinucleotides)
variants = cbind(variants, variantdf)

print(length(colnames(variants)[10:length(colnames(variants))]))
print(length(colnames(variants)[10:length(colnames(variants))]))/length(dinucleotidePropertyNames)

colnames(variants)[10:length(colnames(variants))] = c(paste("1", dinucleotidePropertyNames, sep = "_"), paste("2", dinucleotidePropertyNames, sep = "_"), 
                                         paste("3", dinucleotidePropertyNames, sep = "_"), paste("4", dinucleotidePropertyNames, sep = "_"))

# only keep one column if they are duplicated
variants <- variants[, !duplicated(colnames(variants), fromLast = TRUE)]

variants = variants %>% 
  rename(
    "pos" = "start",
    )

variants = variants[, -3]
print(paste(featureOutputDir,"dinucleotideProperties.txt", sep = ""))
write.table(variants, paste(featureOutputDir,"dinucleotideProperties.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")


