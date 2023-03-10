#%%
import pandas as pd
import os
from liftover import get_lifter
os.chdir("/Users/uw20204/Downloads")
df = pd.read_csv("GRCh37.bed", sep = "\t", header=None)
df
#%%
converter = get_lifter('hg19', 'hg38')
def getLiftover(variant):
    try:
        chrom = df[0][variant]
        pos = df[1][variant]
        refAllele = df[3][variant]
        altAllele = df[4][variant]
        posNew = converter[chrom][pos]
        chrom = posNew[0][0]
        pos1 = int(posNew[0][1])
        pos2 = int(posNew[0][1])
        list = [chrom, pos1, pos2, refAllele, altAllele]
    except:
        list = []
    return(list)

dfLiftover = [getLiftover(variant) for variant in range(0, len(df))]
dfLiftover = pd.DataFrame(dfLiftover)
dfLiftover = dfLiftover[dfLiftover[1].notna()]
dfLiftover = dfLiftover[dfLiftover[2].notna()]
dfLiftover[1] = dfLiftover[1].astype(int)
dfLiftover[2] = dfLiftover[2].astype(int)
dfLiftover = dfLiftover.fillna(0)

dfLiftover.to_csv("/Users/uw20204/Desktop/PhD/GRCh38.bed", sep = "\t", header = None, index=None)
