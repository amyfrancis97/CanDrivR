# Get amino acid substitution matrices
# Extact score for each amino acid change in each matrix

install.packages("Peptides", dependencies=TRUE)
library(devtools)
install_github("dosorio/Peptides")
library(Peptides)
data("AAdata")
install.packages("bios2mds", dependencies=TRUE)
library(bios2mds)
data(sub.mat)

# Upload amino acid changes associated with variants
AA=read.table("/Users/uw20204/Desktop/20221110/AAforRPeptidesPackage.txt", sep = ",", header = TRUE)

# Drop any variants where amino acids are unkown in the dataset
AA = AA[AA['AA1'] != "-", ]
AA = AA[AA['AA2'] != "-", ]

# For each variant, get the score from each amino acid matrix
for(subMatrix in as.list(names(sub.mat))){
  getSubstitMatScores = function(variantRow){
    res = data.frame(sub.mat[subMatrix])[AA[variantRow, 2], gsub(" ", "",paste(names(sub.mat[subMatrix]), '.', AA[variantRow, 3])) ]
    return(res)
  }
  test <- lapply(1:nrow(AA), getSubstitMatScores)
  test[sapply(test, is.null)] <- NA
  test = unlist(test, recursive = TRUE, use.names = TRUE)
  AA = cbind(AA, test)
}

# Write matrix to CSV
colnames(AA) = c(colnames(AA[, 1:5]), names(sub.mat))
AA = AA[-5]
write.csv(AA, "/Users/uw20204/Desktop/20221110/AASubstMatrices.txt", quote = FALSE, row.names = FALSE)
AA