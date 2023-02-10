# Import the AAdata matrix from the peptides package
# The matrix contains different amino acid properties 
# Relating to each of the 20 amino acids
.libPaths("/bp1/mrcieu1/users/uw20204/paper1/features/RpackageLib") 


#install.packages("Peptides", dependencies=TRUE,repos = "http://cran.us.r-project.org")
#install.packages("devtools",repos = "http://cran.us.r-project.org")
library('devtools')
#install_github("dosorio/Peptides", dependencies = TRUE)
#install.packages("bios2mds", dependencies=TRUE,repos = "http://cran.us.r-project.org")
library('bios2mds')
library('Peptides')
data("AAdata")
source("config.R")

# # Upload amino acid changes associated with variants from the VEP AA output
AA=read.table(vepAA , sep = "\t", header = TRUE)

# Drop any variants where amino acids are unkown in the dataset
AA = AA[AA['WT_AA'] != "-", ]
AA = AA[AA['mutant_AA'] != "-", ]

# For each variant, pull out the scores for both the wild type and mutant amino acid
load("AAindex.rda")

getAAExtraProperties = function(variant){
  datalist2=list()
  datalist = c()
  for(i in c("WT_AA", "mutant_AA")){
    if( AA[variant, i] %in% colnames(AAindex)){
      data = AAindex[, AA[variant, i]]
    }else{
      data = rep(NA,nrow(AAindex))
    }
    datalist = c(datalist, data)
    datalist2[[i]] = data
  }
  return(datalist)
  
}

# Carry out function to retrieve amino acid properties for each variant
x = 1:nrow(AA)
lists=lapply(x, getAAExtraProperties)

# Melt lists of variants into a dataframe
df = do.call(rbind, lists)

# Write CSV
# Resultant table are two concatenated vectors
# Vector 1: properties for wild type amino acid
# Vector 2: properties for mutant amino acid
# If synonymous, then the two vectors are identical
colnames(df) = c(paste("WT_AA", AAindex$name, sep = "_"), paste("mutant_AA", AAindex$name, sep = "_"))
df = cbind(AA[, 1:5], df)
name = paste(featureOutputDir,"AAproperties.txt", sep = "")
write.table(df, name, quote = FALSE, row.names = FALSE, sep = "\t")
