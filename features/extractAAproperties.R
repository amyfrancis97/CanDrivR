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
colnames(df) = c(paste("AA1", AAindex$name, sep = "_"), paste("AA2", AAindex$name, sep = "_"))
AAdf2 = cbind(AA[, 1:4], df)
write.table(AAdf2, "/Users/uw20204/Desktop/PhD/VEP_web_aminoAcid_properties2.txt", quote = FALSE, row.names = FALSE, sep = "\t")


