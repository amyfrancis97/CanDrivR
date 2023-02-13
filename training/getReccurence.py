# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np

os.chdir("/Users/uw20204/Downloads")
dfList = []
# reads in dataframe in chunks
with pd.read_csv("CosmicMutantExport.tsv", on_bad_lines='skip', chunksize=100000, sep = "\t", encoding= 'unicode_escape') as reader:
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        chunk = chunk[chunk['HGVSG'].str.contains(">", na=False)] # Filter for substitutions only
        chunk.insert(0, 'chrom', chunk['HGVSG'].str.split(":", expand = True)[0])
        chunk.insert(1, 'start', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(2, 'end', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(3, 'ref_allele', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[-1:])
        chunk.insert(4, 'alt_allele',chunk['HGVSG'].str.split(">", expand = True)[1])
        dfList.append(chunk)
dfList = pd.concat(dfList)
# drop variants whos location is not properly defined
dfList = dfList.drop((dfList[dfList['chrom'].isna()].index), axis = 0)
dfList['R'] = dfList.groupby(['HGVSG'])['ID_sample'].transform('nunique').tolist()
dfList.to_csv("CosmicMutantExportNew.bed", sep = "\t", index = None) # Exports as a bed file
