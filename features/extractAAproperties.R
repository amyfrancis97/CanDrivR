# Import the AAdata matrix from the peptides package
# The matrix contains different amino acid properties 
# Relating to each of the 20 amino acids

install.packages("Peptides", dependencies=TRUE)
install.packages("devtools")
library('devtools')
install_github("dosorio/Peptides", dependencies = TRUE)
install.packages("bios2mds", dependencies=TRUE)
library('bios2mds')
library('Peptides')
data("AAdata")

# Upload amino acid changes associated with variants
AA=read.table("/Users/uw20204/Desktop/20221110/AAforRPeptidesPackage.txt", sep = ",", header = TRUE)

# Drop any variants where amino acids are unkown in the dataset
AA = AA[AA['AA1'] != "-", ]
AA = AA[AA['AA2'] != "-", ]

# For each variant, pull out the scores for both the wild type and mutant amino acid
getAAProperties = function(variant){
  aminoAcidWT = AA[variant, 2]
  aminoAcidMutant = AA[variant, 3]
  # get the value of the element that matches the first amino acid
  WTMproperties = c()
  for(aminoAcidRes in c(aminoAcidWT, aminoAcidMutant)){
    propertiesAAWT = c()
    for(property in 1:length(AAdata)){
      feature = c()
      for(annotation in 1:length(AAdata[[property]])){
        if(aminoAcidRes %in% names(AAdata[[property]][[annotation]])){
          hydrophobicitydf = AAdata[[property]][[annotation]][names(AAdata[[property]][[annotation]]) == aminoAcidRes][[1]]
          feature = c(feature, hydrophobicitydf)
        }else{
          feature = c(feature, "NA")
        }
        
      }
      propertiesAAWT = c(propertiesAAWT, feature)
    }
    var = c(WTMproperties, propertiesAAWT)
  }
  df <- t(data.frame(var))
  return(variantdf)
}

# Carry out function to retrieve amino acid properties for each variant
variantdfapply <- lapply(1:nrow(AA), getAAProperties)

# Melt lists of variants into a dataframe
variantdf = do.call(rbind.data.frame, variantdfapply)

# Write CSV
# Resultant table are two concatenated vectors
# Vector 1: properties for wild type amino acid
# Vector 2: properties for mutant amino acid
# If synonymous, then the two vectors are identical
variantdf = cbind(AA[, 1:4], variantdf[4:ncol(variantdf)])
write.csv(variantdf, "/Users/uw20204/Desktop/20221110/VEP_web_aminoAcid_properties.txt", quote = FALSE, row.names = FALSE)

