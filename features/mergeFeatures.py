import pandas as pd
import ast
import re
import os
from functools import reduce 
import glob
from config import *

dfList = []
for i in glob.glob(featureOutputDir + "*.bedGraph"):
    print(i)
    dfList.append(pd.read_csv(i,  low_memory=False, sep = "\t", header = None, names = ["chrom", "pos", "ref_allele", "alt_allele", "driver_stat", i.split(".")[0] + i.split(".")[1] + ".score"]))

print(dfList[0].head())
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"],
                                            how='outer'), dfList)
print(df_merged)
