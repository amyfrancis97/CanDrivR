#%%
# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np
#%%
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

#%%

os.chdir("/Users/uw20204/Downloads")
dfList = []
# reads in dataframe in chunks
with pd.read_csv("CosmicMutantExportSubs.csv", on_bad_lines='skip', chunksize=300000, sep = ",") as reader:
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        chunk.insert(0, 'chrom', chunk['HGVSG'].str.split(":", expand = True)[0])
        chunk.insert(1, 'start', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(2, 'end', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
        chunk.insert(3, 'ref_allele', chunk['HGVSG'].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[-1:])
        chunk.insert(4, 'alt_allele',chunk['HGVSG'].str.split(">", expand = True)[1])
        dfList.append(chunk)
dfList = pd.concat(dfList)
#%%
os.chdir("/Users/uw20204/Downloads")
dfList = []
# reads in dataframe in chunks
with pd.read_csv("CosmicMutantExportNew.bed", on_bad_lines='skip', chunksize=300000, sep = "\t") as reader:
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        dfList.append(chunk)
dfList = pd.concat(dfList).reset_index(drop = True)
#%%
dfList
#%%
dfList['R'] = dfList.groupby(['HGVSG'])['ID_sample'].transform('nunique').tolist()
#%%
# drop variants whos location is not properly defined
dfList = dfList.drop((dfList[dfList['chrom'].isna()].index), axis = 0)

#%%
print(dfList.head())
#%%
dfList.to_csv("CosmicMutantExportNew.bed", sep = "\t", index = None) # Exports as a bed file
   
# %%
dfList.groupby(['HGVSG','ID_sample']).transform('size')
# %%
len(dfList['ID_sample'].unique())
# %%
dfList.groupby(['HGVSG'])['Primary_histology'].transform('nunique').tolist()
# %%
dfList.columns
# %%
dfList['Gene_CDS_length'].mean()
# %%
