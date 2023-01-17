# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np

os.chdir("/Users/uw20204/Downloads")
dfList = []
# reads in dataframe in chunks
with pd.read_csv("CosmicMutantExportSubs2.csv", on_bad_lines='skip', chunksize=300000, sep = "\t") as reader:
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        chunk = chunk.drop("X", axis =1)
        ref = chunk['POSREF'].str.strip().str[-1].tolist() # Extracts ref allele from REFPOS column
        chunk.insert(2, 'REF', ref)
        chunk['ALT'] = chunk['ALT'].str.strip() # Strips of all white spaces
        chunk['POSREF'] = chunk['POSREF'].str[:-1]
        chunk = chunk.rename(columns = {'POSREF': 'POS'})
        dfList.append(chunk)
dfList = pd.concat(dfList)
dfList = dfList.groupby(dfList.columns.tolist(),as_index=False).size() # Gets the reccurence of each variant
dfList.to_csv("CosmicMutantExport.bed", header = None, sep = "\t", index = None) # Exports as a bed file
