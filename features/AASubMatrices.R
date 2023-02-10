.libPaths("/bp1/mrcieu1/users/uw20204/paper1/features/RpackageLib") 

# Get amino acid substitution matrices
# Extact score for each amino acid change in each matrix
#install.packages("Peptides", dependencies=TRUE,repos = "http://cran.us.r-project.org")
library(devtools)
#install_github("dosorio/Peptides")
library(Peptides)
data("AAdata")
#install.packages("bios2mds", dependencies=TRUE,repos = "http://cran.us.r-project.org")
library(bios2mds)
data(sub.mat)
library("tidyr")
source("config.R")

# # Upload amino acid changes associated with variants from the VEP AA output
AA=read.table(vepAA , sep = "\t", header = TRUE)

# Drop any variants where amino acids are unkown in the dataset
AA = AA[AA['WT_AA'] != "-", ]
AA = AA[AA['mutant_AA'] != "-", ]

# For each variant, get the score from each amino acid matrix
for(subMatrix in as.list(names(sub.mat))){
  getSubstitMatScores = function(variantRow){
    res = data.frame(sub.mat[subMatrix])[AA[variantRow, 'WT_AA'], gsub(" ", "",paste(names(sub.mat[subMatrix]), '.', AA[variantRow, 'mutant_AA'])) ]
    return(res)
  }
  test <- lapply(1:nrow(AA), getSubstitMatScores)
  test[sapply(test, is.null)] <- NA
  test = unlist(test, recursive = TRUE, use.names = TRUE)
  AA = cbind(AA, test)
}

# Write matrix to CSV
colnames(AA) = c(colnames(AA[, 1:7]), names(sub.mat))
AA = AA[ , !(names(AA) %in% c("WT_AA", "mutant_AA"))]
name = paste(featureOutputDir,"AASubstMatrices.txt", sep = "")
write.table(AA, name, quote = FALSE, row.names = FALSE, sep = "\t")

