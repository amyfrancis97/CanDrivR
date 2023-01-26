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
AA=read.table("/Users/uw20204/Desktop/PhD/AAforRPeptidesPackage.txt", sep = ",", header = TRUE)
# Drop any variants where amino acids are unkown in the dataset
AA = AA[AA['AA1'] != "-", ]
AA = AA[AA['AA2'] != "-", ]
################### didnt use this in the end as everything is included in the second table
getAminoAcidProperties = function(AA1, AA2){
  aaPropertyChange = list()
  for(aa in c(AA1, AA2)){
    propertyDfFull2 = c()
    for(property in AAdata){
      propertyDfFull = data.frame()
      for(subproperty in 1:length(property)){
        feature = names(property)[subproperty]
        if(aa %in% names(property[[subproperty]])){
          values = property[[subproperty]][names(property[[subproperty]]) == aa][[1]]
        }else{
          values = NA
        }
        df = data.frame(feature, values)
        propertyDfFull = rbind(propertyDfFull, df)
      }
      propertyDfFull2 = rbind(propertyDfFull2, propertyDfFull)
    }
    aaPropertyChange[[aa]] = propertyDfFull2
  }
  return(aaPropertyChange)
}

AAdf1 = list()
for(row in 1:nrow(AA)){
  AA1 = AA[row, 'AA1']
  AA2 = AA[row, 'AA2']
  aaPropertyChange = getAminoAcidProperties(AA1, AA2)
  if(length(aaPropertyChange) == 2){
    AAdf = merge(aaPropertyChange[[1]], aaPropertyChange[[2]], by = "feature", all = TRUE)
  }else{
    AAdf = merge(aaPropertyChange[[1]], aaPropertyChange[[1]], by = "feature", all = TRUE)
  }
  AAdf1[[row]] = AAdf
}

AAdf2 = data.frame(matrix(ncol = length(c(paste("AA1",AAdf1[[1]][, 'feature'],sep="_"), paste("AA2",AAdf1[[1]][, 'feature'],sep="_")))))
colnames(AAdf2) = c(paste("AA1",AAdf1[[1]][, 'feature'],sep="_"), paste("AA2",AAdf1[[1]][, 'feature'],sep="_"))

for(i in 1:length(AAdf1)){
  AAdf2[i, ] = c(AAdf1[[i]][, 2], AAdf1[[i]][, 3])
}

AAdf2 = cbind(AA[, 1:4], AAdf2)
####################################################
########## GET OTHER PROPERTIES ########
load("/Users/uw20204/Downloads/AAindex.rda")

getAAExtraProperties = function(variant){
  datalist2=list()
  datalist = c()
  for(i in c("AA1", "AA2")){
    if( AAdf2[variant, i] %in% colnames(AAindex)){
      data = AAindex[, AAdf2[variant, i]]
    }else{
      data = rep(NA,nrow(AAindex))
    }
    datalist = c(datalist, data)
    datalist2[[i]] = data
  }
  return(datalist)

}

x = 1:nrow(AAdf2)
lists=lapply(x, getAAExtraProperties)
df = do.call(rbind, lists)
colnames(AAdf2)
colnames(df) = c(paste("AA1", AAindex$name, sep = "_"), paste("AA2", AAindex$name, sep = "_"))
AAdf2 = cbind(AA[, 1:4], df)
AAdf2[1, ]
write.table(AAdf2, "/Users/uw20204/Desktop/PhD/VEP_web_aminoAcid_properties2.txt", quote = FALSE, row.names = FALSE, sep = "\t")


