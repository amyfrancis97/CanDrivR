# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np
from config import *

dfList = []
# reads in dataframe in chunks
with pd.read_csv(cosmicMutantExportCSV, on_bad_lines='skip', chunksize=100000, sep = ",", encoding= 'unicode_escape', header = None) as reader:
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        chunk = chunk[chunk[36].str.contains(">", na=False)] # Filter for substitutions only
        chunk.insert(0, 'chrom', chunk[36].str.split(":", expand = True)[0])
        chunk.insert(1, 'start', chunk[36].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(2, 'end', chunk[36].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(3, 'ref_allele', chunk[36].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[-1:])
        chunk.insert(4, 'alt_allele',chunk[36].str.split(">", expand = True)[1])
        dfList.append(chunk)
dfList = pd.concat(dfList)
print(dfList.head())
# drop variants whos location is not properly defined
dfList = dfList.drop((dfList[dfList['chrom'].isna()].index), axis = 0)
#dfList['R'] = dfList.groupby([36])[36].transform('count')
#dfList['R'] = dfList.groupby([36])[5].transform('nunique').tolist()
#dfList = dfList.drop_duplicates(subset=[35], keep = 'first')
dfList = dfList.groupby(['chrom', 'start', 'end', 'ref_allele', 'alt_allele'], as_index=False).size() # Gets the reccurence of each variant
print(dfList.head())
dfList.to_csv("CosmicMutantExportNew.bed", sep = "\t", index = None) # Exports as a bed file
