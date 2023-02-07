# Get DNA shapes of 10 BP regions overlapping with variant

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("BSgenome")
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
BiocManager::install('DNAshapeR')
install.packages('dplyr')
library('BSgenome')
library('BSgenome.Hsapiens.UCSC.hg38')
library('DNAshapeR')
library('GenomicRanges')
library('dplyr')
library('stringr')
from config import *

# Import variants for shaping
variants=read.table(variants, sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "driver_stat")

# Get the desired base pair range for DNA shape
variants[2] = variants[2]-1
variants[3] = variants[3]+1

# Make a GRRanges object
variants = makeGRangesFromDataFrame(variants)

# Get the 10bp fasta for each variant
getFasta(variants, BSgenome = Hsapiens, width = 3, filename = "/Users/uw20204/Desktop/PhD/VariantDinucleotides.fa")

# Import variants for shaping
variants=read.table("/Users/uw20204/Desktop/PhD/filteredRGreater5.sorted.bed", sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "driver_stat")
VariantDinucleotideWTSeq=read.table("/Users/uw20204/Desktop/PhD/VariantDinucleotides.fa")
toDelete <- seq(1, nrow(VariantDinucleotideWTSeq), 2)
variants = cbind(variants, VariantDinucleotideWTSeq[ -toDelete ,])
VariantDinucleotideWTSeq

getMutantTrinucleotides = function(variantRow){
  mutantTrinucleotides = paste(substr(variants[variantRow, 7], 1, 1), variants[variantRow, 5], substr(variants[variantRow, 7], 3, 3), sep = "")
  return(mutantTrinucleotides)
}

# Carry out function to retrieve mutant trinucleotides for each variant
variantdfapply <- lapply(1:nrow(variants), getMutantTrinucleotides)
variantdfapply
# Melt lists of variants into a dataframe
variantdf = do.call(rbind.data.frame, variantdfapply)
variants = cbind(variants, variantdf)
colnames(variants) = c("chrom", "start", "end", "ref_allele", "alt_allele", "driver_status", "WTtrinuc", "mutTrinuc")

# Read in dinucleotide properties
dinucleotideProperty=read.csv("/Users/uw20204/Desktop/PhD/dinucleotidePropertyTable.csv")

# Get names of dinucleotide properties
dinucleotidePropertyNames = apply(dinucleotideProperty['PropertyName'],2,function(x)gsub('\\s+', '_',x))



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
dinucleotides[2]
# Melt lists of variants into a dataframe
variantdf = do.call(rbind.data.frame, dinucleotides)
variants = cbind(variants, variantdf)

colnames(variants)[9:length(colnames(variants))] = c(paste("1", dinucleotidePropertyNames, sep = "_"), paste("2", dinucleotidePropertyNames, sep = "_"), 
                                         paste("3", dinucleotidePropertyNames, sep = "_"), paste("4", dinucleotidePropertyNames, sep = "_"))


write.csv(variants, featureOutputDir + "dinucleotideProperties.txt", quote = FALSE, row.names = FALSE)

