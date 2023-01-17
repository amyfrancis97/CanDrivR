#%%
import pandas as pd
import ast
import re
import os
#%%
os.chdir('/Users/uw20204/Desktop/20221110')
chunksize = 10000
appended_data = []
for chunk in pd.read_csv("detail_LS_annotationWithPositions.txt", sep = "\t", header = None, chunksize=chunksize):
    LS_annotation = chunk
    LS_annotation['vepID'] = LS_annotation[0] + "_" + LS_annotation[1].astype('string') + "_" + LS_annotation[3] + "/" + LS_annotation[4] 

    SpliceJunctionDf = pd.DataFrame(columns=['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'numSpliceJunctions', 'distancesToLinearSpliceJunctionsSum'] + LS_annotation[13].unique().tolist())
    def getSpliceJunctionInfo(i):
        numSpliceJunctions = len(LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][9].unique())
        uniqueGeneTypes = LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][13].unique().tolist()
        distancesToLinearSpliceJunctionsSum = LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][25].sum()
        SpliceJunctionDf.loc[i, "chrom"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][0].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "vepID"] =  LS_annotation['vepID'].unique().tolist()[i]
        SpliceJunctionDf.loc[i, "pos"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][1].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "ref_allele"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][3].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "alt_allele"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][4].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "driver_status"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][5].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "numSpliceJunctions"] =  numSpliceJunctions
        SpliceJunctionDf.loc[i, "distancesToLinearSpliceJunctionsSum"] =  distancesToLinearSpliceJunctionsSum
        SpliceJunctionDf.loc[i, uniqueGeneTypes] =  1
        return SpliceJunctionDf

    # splits the dataframe up to handle less rows at a time
    input = range(0, len(LS_annotation['vepID'].unique()))
    n = 10000
    output = [input[i:i+n] for i in range(0, len(input), n)]

    for chunk in output:
        [getSpliceJunctionInfo(i) for i in chunk]
    
    # store DataFrame in list
    appended_data.append(SpliceJunctionDf)
appended_data = pd.concat(appended_data)
appended_data1 = appended_data.drop_duplicates(subset='vepID', keep="first")
appended_data1 = appended_data1.fillna(0)
LS_annotation = appended_data1
LS_annotation = LS_annotation.drop("no", axis = 1)
LS_annotation.to_csv("LS_annotation.txt", sep = "\t", index=None)
# %%
os.chdir('/Users/uw20204/Downloads')
pd.read_csv("BS_overlaps.txt", sep = "\t")
#%%
chunksize = 10000
appended_data = []
for chunk in pd.read_csv("BS_overlaps.txt", sep = "\t", header = None, chunksize=chunksize):
    LS_annotation = chunk
    LS_annotation['vepID'] = LS_annotation[0] + "_" + LS_annotation[1].astype('string') + "_" + LS_annotation[3] + "/" + LS_annotation[4] 
    SpliceJunctionDf = pd.DataFrame(columns=['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'numSpliceJunctions', 'distancesToBSJunctionsSum'] + LS_annotation[13].unique().tolist())
    def getSpliceJunctionInfo(i):
        numSpliceJunctions = len(LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][9].unique())
        uniqueGeneTypes = LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][13].unique().tolist()
        distancesToLinearSpliceJunctionsSum = LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][18].sum()
        SpliceJunctionDf.loc[i, "chrom"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][0].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "vepID"] =  LS_annotation['vepID'].unique().tolist()[i]
        SpliceJunctionDf.loc[i, "pos"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][1].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "ref_allele"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][3].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "alt_allele"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][4].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "driver_status"] =  LS_annotation[LS_annotation['vepID'] == LS_annotation['vepID'].unique().tolist()[i]][5].reset_index(drop=True)[0]
        SpliceJunctionDf.loc[i, "numSpliceJunctions"] =  numSpliceJunctions
        SpliceJunctionDf.loc[i, "distancesToBSJunctionsSum"] =  distancesToLinearSpliceJunctionsSum
        SpliceJunctionDf.loc[i, uniqueGeneTypes] =  1
        return SpliceJunctionDf
    
    # splits the dataframe up to handle less rows at a time
    input = range(0, len(LS_annotation['vepID'].unique()))
    n = 10000
    output = [input[i:i+n] for i in range(0, len(input), n)]

    for chunk in output:
        [getSpliceJunctionInfo(i) for i in chunk]
    
    # store DataFrame in list
    appended_data.append(SpliceJunctionDf)
# %%
appended_data = pd.concat(appended_data)
# %%
appended_data1 = appended_data.drop_duplicates(subset='vepID', keep="first")
appended_data1 = appended_data1.fillna(0)
LS_annotation = appended_data1
#%%
LS_annotation = LS_annotation.drop(".", axis = 1)
LS_annotation = LS_annotation.drop("nan", axis = 1)
#%%
LS_annotation = LS_annotation.rename(columns = {"distancesToLinearSpliceJunctionsSum": "distancesToBSJunctionsSum"})

#%%
LS_annotation.to_csv("/Users/uw20204/Desktop/20221110/BS_annotation.txt", sep = "\t", index=None)
# %%
LS_annotation
# %%
