# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np
from config import *

dfList = []
# reads in dataframe in chunks
with pd.read_csv(gnomadFile, on_bad_lines='skip', chunksize=50000, sep = "\t", header = None) as reader:
    for chunk in reader:
        df = chunk.reset_index(drop=True)
        df = df.loc[(df[5].str.len() == 1) & (df[4].str.len() == 1)]
        try:
            df = df.loc[(df[8].str.split("non_cancer_AF=", expand = True)[1].str.split(";", expand = True)[0].astype("float") > 0.01) |
            (df[8].str.split("controls_AF_raw=", expand = True)[1].str.split(";", expand = True)[0].astype("float") > 0.01)]
        except:
            df = df.loc[df[8].str.split("controls_AF_raw=", expand = True)[1].str.split(";", expand = True)[0].astype("float") > 0.01]
        dfList.append(df)
dfList = pd.concat(dfList)
dfList.to_csv(outputDir + "gnomadGreater0.01AF.bed", header = None, index = None, sep = "\t")
