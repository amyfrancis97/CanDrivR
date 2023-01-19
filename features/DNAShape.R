# Get DNA shapes of 10 BP regions overlapping with variant
# Gets the shape of the wild-type regions with the ref allele

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("BSgenome")
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
BiocManager::install('DNAshapeR')
library('BSgenome')
library('BSgenome.Hsapiens.UCSC.hg38')
library('DNAshapeR')
library('GenomicRanges')

# Import variants for shaping
variants=read.table("/Users/uw20204/Desktop/PhD/filteredRGreater5.sorted.bed", sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "driver_stat")
variants = variants[1:5, ]
# Get the desired base pair range for DNA shape
variants[2] = variants[2]-5
variants[3] = variants[3]+5
variants
# Make a GRRanges object
variants = makeGRangesFromDataFrame(variants)

# Get the 10bp fasta for each variant
getFasta(variants, BSgenome = Hsapiens, width = 10, filename = "/Users/uw20204/Desktop/PhD/hg38.fa")
fn <- "/Users/uw20204/Desktop/PhD/hg38.fa"
fn
# Get the shape properties of each position for each variant
pred <- getShape(fn)
names(pred)
install.packages('tidyverse')
library(tidyverse)

# Reduce all of the properties into a single matrix
dnaShape = Reduce("cbind", pred)
variants=read.table("/Users/uw20204/Desktop/PhD/filteredRGreater5.sorted.bed", sep = "\t")
colnames(variants) =  c("chrom", "start", "end", "ref", "alt", "driver_stat")
dnaShape = cbind(variants, dnaShape)

colnames(dnaShape) = c(colnames(dnaShape)[1:6], paste(1:10, "MGW", sep = "_"), paste(1:9, "HelT", sep = "_"), paste(1:10, "ProT", sep = "_"),
                 paste(1:9, "Roll", sep = "_"), paste(1:10, "EP", sep = "_"))


# Write DNA shape properties to CSV
write.csv(dnaShape, "/Users/uw20204/Desktop/PhD/dnaShape.txt", quote = FALSE, row.names = FALSE)

