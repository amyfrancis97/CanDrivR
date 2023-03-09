# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np
from config import *
import sys 

if __name__ == "__main__":
    dfList = []
    col = int(sys.argv[2])-1
    file = sys.argv[1]
    # reads in dataframe in chunks
    with pd.read_csv(file + ".csv", on_bad_lines='skip', chunksize=100000, sep = ",", encoding= 'unicode_escape', header = None) as reader:
        for chunk in reader:
            chunk = chunk.reset_index(drop=True)
            chunk = chunk[chunk[col].str.contains(">", na=False)] # Filter for substitutions only
            chunk.insert(0, 'chrom', chunk[col].str.split(":", expand = True)[0])
            chunk.insert(1, 'start', chunk[col].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
            chunk.insert(2, 'end', chunk[col].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[:-1])
            chunk.insert(3, 'ref_allele', chunk[col].str.split(":", expand = True)[1].str.split("g.", expand = True)[1].str.split(">", expand = True)[0].str[-1:])
            chunk.insert(4, 'alt_allele',chunk[col].str.split(">", expand = True)[1])
            dfList.append(chunk)
    dfList = pd.concat(dfList)
    print(dfList.head())
    # drop variants whos location is not properly defined
    dfList = dfList.drop((dfList[dfList['chrom'].isna()].index), axis = 0)
    dfList = dfList.groupby(['chrom', 'start', 'end', 'ref_allele', 'alt_allele'], as_index=False).size() # Gets the reccurence of each variant
    print(dfList.head())
    dfList.to_csv(sys.argv[1] + ".bed", sep = "\t", index = None) # Exports as a bed file
