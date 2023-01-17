# Get spice data 
wget http://www.rjunbase.org/rJunBase/download/downloadData?fileName=detail_FS_annotation.txt

# Reformat data
python X.py

# Get closest feature to each 
bedtools closest -d  -a /Users/uw20204/Desktop/20221110/filteredRGreater5.sorted.txt -b BS_annotation.txt > BS_overlaps.txt

#######
# Get alternative splice junctions 
wget http://www.rjunbase.org/rJunBase/download/downloadData?fileName=Alternative_Splice_Junctions.tar.gz

tar -xf Alternative_Splice_Junctions.tar.gz

python X.py

# Get closest feature to each 
bedtools closest -d  -a /Users/uw20204/Desktop/20221110/filteredRGreater5.sorted.txt -b alternative_splicing_annotation.txt > AS_overlaps.txt

python Y.py


