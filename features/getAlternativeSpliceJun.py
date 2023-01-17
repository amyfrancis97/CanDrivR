#%%
import pandas as pd
import ast
import re
import os
#%%
os.chdir('/Users/uw20204/Downloads')
chunksize = 10000
appended_data = []
for chunk in pd.read_csv("Alternative_Splice_Junctions.csv", sep = ",", chunksize=chunksize):
    chunk = chunk.reset_index(drop = True)
    chrom = chunk['Junction location'].str.split(":",expand = True)[0].str.strip()
    pos1 = chunk['Junction location'].str.split(":",expand = True)[1].str.split("|",expand = True)[0].str.strip()
    pos2 = chunk['Junction location'].str.split(":",expand = True)[1].str.split("|",expand = True)[1].str.strip()
    chunk.insert(0, 'chrom', chrom)
    chunk.insert(1, 'pos1', pos1)
    chunk.insert(2, 'pos2', pos2)

    chunk1 = chunk['Alternative splicing types'].str.split(";", expand = True)
    altSplicetypes = []
    for col in chunk1.columns:
        altSplicetypes.append(chunk1[col].str.split(":", expand = True)[0])
    altSplicetypes = pd.concat(altSplicetypes, axis = 1)
    def getAltSplice(variant):
        altSpliceList = altSplicetypes.iloc[variant, :].unique().tolist()
        res = [i for i in altSpliceList if i is not None]
        chunk.loc[variant, res] = 1
    [getAltSplice(variant) for variant in range(0, len(chunk))]
    chunk = chunk.drop(["Junction location", "Alternative splicing types"], axis =1)
    # store DataFrame in list
    appended_data.append(chunk)

appended_data = pd.concat(appended_data)
indicesToChange = appended_data[appended_data['pos2'] < appended_data['pos1']].index.tolist()
appended_data.loc[indicesToChange, 'pos1']= appended_data.loc[indicesToChange, 'pos2']
appended_data.loc[indicesToChange, 'pos2']= appended_data.loc[indicesToChange, 'pos1']
appended_data['pos1'] = appended_data['pos1'].astype(int)
appended_data['pos2'] = appended_data['pos2'].astype(int)
appended_data['chrom'] = appended_data['chrom'].str.strip()
appended_data = appended_data.sort_values(['chrom', 'pos1', 'pos2'],
              ascending = [True, True, True])
appended_data = appended_data.fillna(0)
appended_data.to_csv("alternative_splicing_annotation.txt", sep = "\t", index=None, header = None)
#%%

# %%
os.chdir('/Users/uw20204/Downloads')
AS_overlaps = pd.read_csv("AS_overlaps.txt", sep = "\t", header = None)
AS_overlaps['vepID'] = AS_overlaps[0] + "_" + AS_overlaps[1].astype('string') + "_" + AS_overlaps[3] + "/" + AS_overlaps[4] 
AS_overlaps = AS_overlaps.drop([6, 7, 8], axis = 1)
altSpliceDf = pd.DataFrame(columns = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'vepID', 'mutex_exons','alt_3prime', 'alt_5prime', 'exon_skip', 'intron_retention'])

def getFinalAltSpliceDf(variant):
    chunked = AS_overlaps[AS_overlaps['vepID'] == AS_overlaps['vepID'].unique()[variant]].reset_index(drop = True)
    chrom = chunked[0][0]
    pos = chunked[1][0]
    refAllele = chunked[3][0]
    altAllele = chunked[4][0]
    driverStatus = chunked[5][0]
    vepID = chunked['vepID'][0]

    if 1.0 in chunked.iloc[:, 6].tolist():
        mutex_exons = 1
    else:
        mutex_exons = 0
    if 1.0 in chunked.iloc[:, 7].tolist():
        alt_3prime = 1
    else:
        alt_3prime = 0
    if 1.0 in chunked.iloc[:, 8].tolist():
        alt_5prime = 1
    else:
        alt_5prime = 0
    if 1.0 in chunked.iloc[:, 9].tolist():
        exon_skip = 1
    else:
        exon_skip = 0
    if 1.0 in chunked.iloc[:, 10].tolist():
        intron_retention = 1
    else:
        intron_retention = 0

    altSpliceDf.loc[variant, :] = [chrom, pos, refAllele, altAllele, driverStatus, vepID, mutex_exons, alt_3prime, alt_5prime, exon_skip, intron_retention]
    return altSpliceDf

dfList = [getFinalAltSpliceDf(variant) for variant in range(0, len(AS_overlaps['vepID'].unique()))]
dfList = pd.concat(dfList)
# %%
appended_data.to_csv("/Users/uw20204/Desktop/20221110/AS_annotation.txt", sep = "\t", index=None)

# %%
appended_data
# %%
