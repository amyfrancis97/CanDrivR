# Get the chromosome and position of the LS annotation file 
# Reformatting so that we can use closest bedtools function in command line
library(stringr)
LS_anno = read.csv("/Users/uw20204/Downloads/detail_LS_annotation.txt", sep = "\t")

getChromPos = function(row){
  chrom = str_split_i(LS_anno[row, 'Junction.location'], ":", 1)
  pos = str_split_i(LS_anno[row, 'Junction.location'], ":", 2)
  pos1 = str_split_i(pos, "\\|", 1)
  pos2 = str_split_i(pos, "\\|", 2)
  return(data.frame(chrom, pos1, pos2))
}

chromPos <- lapply(1:nrow(LS_anno), getChromPos)
df1 = do.call(rbind, chromPos)
LS_anno = cbind(df1, LS_anno)
write.table(LS_anno, "/Users/uw20204/Downloads/detail_LS_annotationWithPositions.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
