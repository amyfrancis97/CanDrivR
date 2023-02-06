#%%
import pandas as pd
import ast
import re
import os

os.chdir('/Users/uw20204/Desktop/PhD/')

gnomadVEP1 = pd.read_csv("gnomadVEP.txt", sep = "\t", low_memory=False)
cosmicVEP1 = pd.read_csv("cosmicVEP.txt", sep = "\t", low_memory=False)

def readCons(cosmicFileName, gnomadFileName):
    cosmic = pd.read_csv(cosmicFileName + ".bed", sep = "\t", header=None, names = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', cosmicFileName.replace("cosmic", "") + "Score"])
    gnomad = pd.read_csv(gnomadFileName + ".bed", sep = "\t", header=None, names = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', gnomadFileName.replace("gnomad", "") + "Score"])
    cosmic['driver_status'] = 1
    gnomad['driver_status'] = 0
    # make column to match VEP output to merge variant scores
    cosmic['vepID'] = cosmic.chrom + "_" + cosmic.pos.astype(str) + "_" + cosmic.ref_allele + "/" + cosmic.alt_allele
    cosmic['vepID'] = cosmic['vepID'].astype('string')
    gnomad['vepID'] = gnomad.chrom + "_" + gnomad.pos.astype(str) + "_" + gnomad.ref_allele + "/" + gnomad.alt_allele
    gnomad['vepID'] = gnomad['vepID'].astype('string')
    df = pd.concat([cosmic, gnomad])
    return df

os.chdir('/Users/uw20204/Desktop/PhD/conservation')
# reading in phyloP scores
phyloP100 = readCons("cosmicPhyloP100", "gnomadPhyloP100")
phyloP30 = readCons("cosmicPhyloP30", "gnomadPhyloP30")
phyloP20 = readCons("cosmicPhyloP20", "gnomadPhyloP20")
phyloP17 = readCons("cosmicPhyloP17", "gnomadPhyloP17")
phyloP4 = readCons("cosmicPhyloP4", "gnomadPhyloP4")
phyloP7 = readCons("cosmicPhyloP7", "gnomadPhyloP7")

# reading in phastcons scores
phastCons100 = readCons("cosmicPhastCons100", "gnomadPhastCons100")
phastCons30 = readCons("cosmicPhastCons30", "gnomadPhastCons30")
phastCons20 = readCons("cosmicPhastCons20", "gnomadPhastCons20")
phastCons17 = readCons("cosmicPhastCons17", "gnomadPhastCons17")
phastCons4 = readCons("cosmicPhastCons4", "gnomadPhastCons4")
phastCons7 = readCons("cosmicPhastCons7", "gnomadPhastCons7")

os.chdir('/Users/uw20204/Desktop/PhD')

# Merging all conservation datasets
cons = phastCons20.merge(phyloP20, on = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', 'driver_status', 'vepID'])
for consFeature in [phyloP100, phastCons100, phastCons4, phyloP4, phastCons7, phyloP7, phastCons30, phyloP30, phyloP17, phastCons17]:
    cons = cons.merge(consFeature, on = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', 'driver_status', 'vepID'])
#%%
consCols = cons.columns.tolist()[8:] + [cons.columns.tolist()[5]]
#%%
# Read in GC content features
dfGC = []
for i in [10000, 2000, 1000, 500, 200, 100]:
    dfGC.append(pd.read_csv(str(i) + "_GC.txt", sep = ",").iloc[:, range(0,7)])


indices = dfGC[0].iloc[:, range(0, 5)]
dfGC = pd.concat(dfGC, axis=1)
dfGC['driver_status'] = 0
dfGC = pd.concat([indices, dfGC.drop(indices.columns, axis = 1)], axis = 1)

dfGC = dfGC.rename(columns = {'reference_allele': 'ref_allele', 'alternate_allele': 'alt_allele'})
dfGC['chrom'] = "chr" + dfGC['chrom'].astype('string')
#%%
GCCols = dfGC.columns.tolist()[5:]
#%%

############# kernel ###########
kernel_5_1 = pd.read_table("5_1_kernel.txt", sep = ",", index_col=0)
kernel_5_2 = pd.read_table("5_2_kernel.txt", sep = ",", index_col=0)
kernel_5_3 = pd.read_table("5_3_kernel.txt", sep = ",", index_col=0)
kernel_4_1 = pd.read_table("4_1_kernel.txt", sep = ",", index_col=0)
kernel_4_2 = pd.read_table("4_2_kernel.txt", sep = ",", index_col=0)
kernel_4_3 = pd.read_table("4_3_kernel.txt", sep = ",", index_col=0)
kernel_3_1 = pd.read_table("3_1_kernel.txt", sep = ",", index_col=0)
kernel_3_2 = pd.read_table("3_2_kernel.txt", sep = ",", index_col=0)
kernel_3_3 = pd.read_table("3_3_kernel.txt", sep = ",", index_col=0)
kernel_2_2 = pd.read_table("2_2_kernel.txt", sep = ",", index_col=0)
kernel_2_1 = pd.read_table("2_1_kernel.txt", sep = ",", index_col=0)
kernel_1_1 = pd.read_table("1_1_kernel.txt", sep = ",", index_col=0)

def changeColNames(kernelDf, kernelFeatureNum):
    kernelDf = kernelDf.rename(columns = {'0':str(kernelFeatureNum) + '_0', '1':str(kernelFeatureNum) + '_1',
    '2':str(kernelFeatureNum) + '_2' , '3':str(kernelFeatureNum) + '_3'})
    return kernelDf

dfKernel = []
kernelFeatures = [kernel_5_3, kernel_5_2, kernel_5_1,
kernel_4_3, kernel_4_2, kernel_4_1, kernel_3_3, kernel_3_2, kernel_3_1,
kernel_2_2, kernel_2_1, kernel_1_1]
for kernelFeatureNum in range(0,len(kernelFeatures)):
    newDf = changeColNames(kernelFeatures[kernelFeatureNum], kernelFeatureNum)
    dfKernel.append(newDf)

dfKernel = pd.concat(dfKernel, axis =1)
dfKernel = dfKernel.loc[:, ~dfKernel.columns.duplicated(keep='first')]

dfKernel['vepID'] = "chr" + dfKernel['chrom'].astype('string') + "_" + dfKernel['pos'].astype('string') + "_" + dfKernel['reference_allele'] + "/" + dfKernel['alternate_allele']
dfKernel = dfKernel.drop(dfKernel.columns.to_list()[0:5], axis = 1)
#%%
kernelCols = dfKernel.columns.tolist()[1:len(dfKernel.columns.tolist())-1]
#%%
AAproperties = pd.read_csv("VEP_web_aminoAcid_properties2.txt", sep = "\t")

#%%

#%%
from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
# define example
data1 = pd.DataFrame(AAproperties['AA1'])
from sklearn.preprocessing import OneHotEncoder
ohe1 = OneHotEncoder(sparse=False)
result1 = ohe1.fit_transform(data1)

# define example
data2 = pd.DataFrame(AAproperties['AA2'])
from sklearn.preprocessing import OneHotEncoder
ohe2 = OneHotEncoder(sparse=False)
result2 = ohe2.fit_transform(data2)

vepAA_autoencoded = pd.concat([pd.DataFrame(result1, columns=pd.DataFrame(["WTAA_" + x for x in ohe1.categories_]).iloc[0, :].tolist()), pd.DataFrame(result2, columns=pd.DataFrame(["mutantAA_" + x for x in ohe2.categories_]).iloc[0, :].tolist())], axis =1)

vepAA_autoencoded.insert(0, "vepID", AAproperties['vepID'])
#%%
AAautoencodedCols = vepAA_autoencoded.columns.tolist()[1:]
#%%
AAproperties = AAproperties.drop(['AA1', 'AA2'], axis = 1)
#%%
AApropertyCols = AAproperties.columns.tolist()[2:]
#%%
ATAC = pd.read_csv("intersectAbsSumit.bed", sep = "\t", header=None)
ATAC = ATAC.drop([2, 6, 7, 8], axis = 1)
ATAC = ATAC.rename(columns = {0: 'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5: 'driver_status', 9:'mean_pileup', 10:'median_pileup', 11:'mean_fold_enrichment', 12: 'median_fold_enrichment', 13: 'distance_feature_to_variant'})
ATAC['vepID'] = ATAC['chrom'] + "_" + ATAC['pos'].astype('string') + "_" + ATAC['ref_allele'] + "/" + ATAC['alt_allele'] 
# takes the mean of any duplicates
ATAC2 = ATAC.groupby('vepID').mean().reset_index()

ATAC2 = ATAC2.rename(columns = {'mean_pileup': 'ATAC_mean_pileup', 'median_pileup': 'ATAC_median_pileup', 'mean_fold_enrichment': 'ATAC_mean_fold_enrichment', 
'median_fold_enrichment':'ATAC_median_fold_enrichment', 'distance_feature_to_variant':'ATAC_distance_feature_to_variant'})
#%%
ATACcols = ATAC2.columns.tolist()[3:]
#%%
TSSdistance = pd.read_csv("intersectTSS.bed", sep = "\t", header = None)
TSSdistance = TSSdistance.drop([2, 6, 7, 8, 9, 10, 11, 12, 13, 14], axis =1)

TSSdistance = TSSdistance.rename(columns = {0: 'chrom', 1: 'pos', 3:'ref_allele', 4: 'alt_allele', 5: 'driver_status', 15: 'distanceTSS'})
#%%
tssDistanceCols = TSSdistance.columns.tolist()[5]
#%%
def uniquenessDf(x, y, z):
    uniqueness = pd.read_csv(x+y, sep = "\t", names = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', (y + 'Score')])
    uniqueness['driver_status'] = z
    return uniqueness

dfUniqueGnomad = pd.DataFrame(columns={'chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance'})
for i in ['k100.Bismap.MultiTrackMappability.bed', 'k100.Umap.Multi.bed']:
    uniquenessGnomad = uniquenessDf('gnomad', i, 0)
    dfUniqueGnomad = dfUniqueGnomad.merge(uniquenessGnomad, how = 'outer')

dfUniqueCosmic = pd.DataFrame(columns={'chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance'})
for i in ['k100.Bismap.MultiTrackMappability.bed', 'k100.Umap.Multi.bed']:
    uniquenessCosmic = uniquenessDf('cosmic', i, 1)
    dfUniqueCosmic = dfUniqueCosmic.merge(uniquenessCosmic, how = 'outer')

dfUnique = pd.concat([dfUniqueCosmic, dfUniqueGnomad], axis = 0)
#%%
uniquenessCols = [dfUnique.columns.tolist()[5]] + [dfUnique.columns.tolist()[7] ]
#%%
vepAA = pd.read_csv("featuresTraining/vepAA.txt", sep = "\t")
vepAA = vepAA.rename(columns = {'AA1': 'WT_AA', 'AA2': 'mutant_AA'})
#%%
vep = pd.read_csv("featuresTraining/vep.txt", sep = "\t")
#cons = pd.read_csv("featuresTraining/cons.txt", sep = "\t")
#%%
vepConseqCols = vep.columns.tolist()[1:len(vep.columns.tolist())-1]
#%%
AASubstMatrix = pd.read_csv("AASubstMatrices.txt", sep = ",")
AASubstMatrix = AASubstMatrix.drop(['AA1', 'AA2'], axis =1)
#%%
AASubstCols = AASubstMatrix.columns.tolist()[2:]
#%%
dnaShape = pd.read_csv("dnaShape.txt", sep = ",")
dnaShape = dnaShape.dropna(axis=1, how='all')
dnaShape = dnaShape.drop(['end', 'driver_stat'], axis = 1)
dnaShape = dnaShape.rename(columns = {'start': 'pos', 'ref' : 'ref_allele', 'alt': 'alt_allele'})
dnaShapeCols = dnaShape.columns.tolist()[4:len(dnaShape.columns.tolist())-1]
L1 = [str(col) + '_DNAshape' for col in dnaShape.columns[4:]]
dnaShape = dnaShape.set_axis(dnaShape.columns.tolist()[:4] + L1, axis=1, inplace=False)
dnaShapeCols = dnaShape.columns.tolist()[4:len(dnaShape.columns.tolist())-1]
dnaShapeCols
#%%
LS_annotation = pd.read_csv("LS_annotation.txt", sep = "\t")

#%%
LSannoCols = LS_annotation.columns.tolist()[6:]
#%%
#AS_annotation = pd.read_csv("AS_annotation.txt", sep = "\t")
#AS_annotation = AS_annotation.drop("pos1")
#AS_annotation = AS_annotation.rename(columns = {'pos2': 'pos'})
#%%
dinucleotideProperties = pd.read_csv("dinucleotideProperties.txt", sep = ",")
dinucleotideProperties = dinucleotideProperties.rename(columns = {'start': 'pos'}).drop(['end', 'WTtrinuc', 'mutTrinuc'], axis =1)
#%%
dinucleotideProperties['driver_status']
#%%
dinucPropCols = dinucleotideProperties.columns.tolist()[5:len(dinucleotideProperties.columns.tolist())-1]

#%%
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + dinucPropCols]

df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + LSannoCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + dnaShapeCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + AASubstCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + vepConseqCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + uniquenessCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + tssDistanceCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + ATACcols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + AApropertyCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + AAautoencodedCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + kernelCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + GCCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + consCols]
#%%
# Merge all of the current feature datasets together
df = cons.merge(vep)
df = df.merge(vepAA_autoencoded, on = 'vepID')
df = df.merge(dfKernel)
df = df.merge(dfGC)
df = df.merge(dinucleotideProperties)
df = df.merge(AAproperties) #, on = ['vepID', 'driver_status']
df = df.merge(AASubstMatrix)
df = df.merge(dnaShape)
df = df.merge(ATAC2)
df = df.merge(TSSdistance)
df = df.merge(dfUnique)
df = df.merge(LS_annotation)

#%%
# drop variants from cosmic dataset that exist in gnomad as we can asssume that these are neutral
duplicates = list(df[df.duplicated('vepID')].sort_values(['chrom', 'pos'])['vepID'])
df = df.drop(df[df.duplicated('vepID', keep = False)][df[df.duplicated('vepID', keep = False)]['driver_status'] == 1].index)
df = df.reset_index(drop = True)

df = df[(df.chrom != 'chrY') & (df.chrom != 'chrX')].reset_index(drop=True)
# Restrict training data to driver that occur in >5 instances/samples
# this will also restrict the negative training data only to variants that overlap with new data points
# 1000 base pair window
df = df[df.reccurance > 5].reset_index(drop=True)

# grouping by groups of chromosomes rather than individual chromosomes
# increases size of test dataset and speeds up the LOCO process
df['grouping'] = 0
for i in ['chr1', 'chr2', 'chr3']:
    df.loc[df['chrom'] == i, 'grouping'] = 1
for i in ['chr4', 'chr5', 'chr6']:
    df.loc[df['chrom'] == i, 'grouping'] = 2
for i in ['chr7', 'chr8', 'chr9']:
    df.loc[df['chrom'] == i, 'grouping'] = 3
for i in ['chr10', 'chr11', 'chr12']:
    df.loc[df['chrom'] == i, 'grouping'] = 4
for i in ['chr13', 'chr14', 'chr15']:
    df.loc[df['chrom'] == i, 'grouping'] = 5
for i in ['chr16', 'chr17', 'chr18']:
    df.loc[df['chrom'] == i, 'grouping'] = 6
for i in ['chr19', 'chr20', 'chr21', 'chr22']:
    df.loc[df['chrom'] == i, 'grouping'] = 7
#%%
df
#%%
df.to_csv("featuresAll.txt", sep = "\t", index=None)
#%%
df
# %%
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance'] + dinucPropCols]

df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + LSannoCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + dnaShapeCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + AASubstCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + vepConseqCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + uniquenessCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + tssDistanceCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + ATACcols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + AApropertyCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + AAautoencodedCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + kernelCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + GCCols]
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + consCols]
##################
#%%
df
#%%
df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + dnaShapeCols]
#%%
res = testLOCO(df[['vepID', 'chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_status', 'reccurance', 'grouping'] + dinucPropCols])
#%%
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

#%%
featureGroups = [dinucPropCols, dnaShapeCols, AASubstCols, LSannoCols, vepConseqCols, uniquenessCols, tssDistanceCols,
ATACcols, AApropertyCols, AAautoencodedCols, kernelCols, GCCols, consCols]
featureScoresres = []
def getScoresPerFeature(featureGroup):
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj][0]

    if isinstance(featureGroup, list):
        featureGroup2 = featureGroup
        name = namestr(featureGroup, globals())
    else:
        name = featureGroup
        featureGroup2 = [featureGroup]
    
    singleDf = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"] + featureGroup2]
    return [name, testLOCO(singleDf)[0], testLOCO(singleDf)[1]]

#%%
featureScores = [getScoresPerFeature(feature) for feature in featureGroups]
featureScores = pd.DataFrame(featureScores)
#%%

featureScores.sort_values(1, ascending=False)[0].tolist()
#%%
singleDf = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", 
"grouping"] + AApropertyCols + AAautoencodedCols + consCols + AASubstCols + dinucPropCols + 
vepConseqCols + dnaShapeCols + kernelCols + LSannoCols + uniquenessCols]
testLOCO(singleDf)
#%%
featureGroups
#%%
a_list= []
resAll = []
resAllStd = []
names = []
for feature in [AApropertyCols, AAautoencodedCols, consCols, AASubstCols, dinucPropCols,
vepConseqCols, dnaShapeCols, kernelCols, LSannoCols, uniquenessCols, GCCols, tssDistanceCols, ATACcols]:
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj][0]

    if isinstance(feature, list):
        a_list = a_list + feature   
        name = namestr(feature, globals())
    else:
        a_list = a_list + [feature]
        name = feature
    names.append(name)
    print(names)
    dfTopFeatures = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"] + a_list]
    res = testLOCO(dfTopFeatures)
    res1 = res[0]
    res2 = res[1]
    resAll.append(res1)
    resAllStd.append(res2)
#%%

import xgboost as xgb
from xgboost import XGBClassifier
from xgboost import plot_importance
from xgboost.sklearn import XGBClassifier

# sklearn packages
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn import svm
from sklearn.utils import shuffle
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn import svm, datasets
from sklearn import metrics
from sklearn.feature_selection import SelectFromModel
from sklearn import metrics   #Additional scklearn functions
from sklearn.model_selection import LeaveOneGroupOut
from numpy import asarray
from sklearn.preprocessing import StandardScaler

# other packages
from ctypes import Structure
import numpy as np
from numpy import mean
from numpy import std
from numpy import loadtxt
import pandas as pd
import sys
import os
from matplotlib import pyplot
from matplotlib.pylab import rcParams
from xml.sax.handler import feature_namespace_prefixes
# plot feature importance using built-in function
from numpy import loadtxt
from xgboost import XGBClassifier
from xgboost import plot_importance
from matplotlib import pyplot
from numpy import loadtxt
from xgboost import XGBClassifier
from xgboost import plot_importance
from matplotlib import pyplot
from numpy import loadtxt
from numpy import sort
from xgboost import XGBClassifier

from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from sklearn import datasets
from sklearn.neighbors import KNeighborsClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split
from sklearn.calibration import calibration_curve
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.svm import LinearSVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_classification
from sklearn.metrics import roc_auc_score
from sklearn.calibration import CalibrationDisplay
import gc
import glob
import os
import random
import time
from datetime import date, datetime

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn import model_selection
from sklearn.metrics import accuracy_score, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

%matplotlib inline
# samples positive dataset to match number of negatives
# may have to repeat this step to train on all positive data?
import warnings
warnings.filterwarnings('ignore') # setting ignore as a parameter

def testLOCO(dataframe):
    mean1 = []
    std1 = []
    for i in range(1, 3):
        # randomly sample dataset for balances classes
        df = pd.concat([dataframe[dataframe.driver_status == 1].sample(len(dataframe[dataframe.driver_status == 0])), dataframe[dataframe.driver_status == 0]]).reset_index(drop=True)
        X = df.drop(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"], axis=1)
        #X = dataframe.drop(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "grouping"], axis=1)
        y = df["driver_status"]
        groups = df["grouping"]
        logo = LeaveOneGroupOut()
        predProbDf1 = []
        indices = []
        scoresXGBMean = []
        scoresXGBStd = []
        scoresXGB = []
        roc_auc1 = []
        test_values_per_fold = []
        SHAP_values_per_fold = [] 
        featureImportance=[]
        keys = []
        rocAll = []
        for train_index, test_index in logo.split(X, y, groups):
            indices.append([list(train_index), (list(test_index))])
            X_train, X_test = X.iloc[train_index, :].reset_index(drop = True), X.iloc[test_index, :].reset_index(drop = True)
            y_train, y_test = y[train_index].reset_index(drop = True), y[test_index].reset_index(drop = True)

            from numpy import asarray
            from sklearn.preprocessing import StandardScaler
            # define standard scaler
            scaler = StandardScaler()
            # transform data
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.fit_transform(X_test)
            xgb_model =  XGBClassifier()
            #xgb_model =  SVC(probability=True)
            # Fit the model with training data and target values
            xgb_model.fit(X_train, y_train)
            # Explain model predictions using shap library:
            #explainer = shap.TreeExplainer(xgb_model)
            #shap_values = explainer.shap_values(X_test)
            #for SHAPs in shap_values:
                #SHAP_values_per_fold.append(SHAPs)
            #for test in X_test:
                #test_values_per_fold.append(test)
            feature_names = X.columns
            featureImportance.append(xgb_model.feature_importances_)
            # Plot summary_plot
            #shap.summary_plot(shap_values, X_test, feature_names = X.columns)
    #####################
            y_pred = xgb_model.predict(X_test)


            predProb = xgb_model.predict_proba(X_test)

            predProbDf = pd.DataFrame(predProb)
            predProbDf.insert(0, 'PredictedLabel', y_pred)
            predProbDf.insert(0, 'ActualLabel', y_test.values)
            predProbDf = predProbDf.rename(columns = {0: 'prob(0)', 1: 'prob(1)'})

            #platt = PlattCalibrator(log_odds=True).fit()
            #platt.fit(predProb[:, 1], y_test.values)
            #platt_probs = platt.predict(xgb_pred_test)
            #platt_probs
            #iso_reg = IsotonicRegression(y_min = 0, y_max = 1, out_of_bounds = 'clip').fit(prob_val, y_val)
            #prob_calibrated = iso_reg.predict(model.predict_proba(X_test)[:, 1])

            ##############
            # Get ROC
            preds = predProb[:,1]
            fpr, tpr, threshold = metrics.roc_curve(y_test, preds)

            roc_auc = metrics.auc(fpr, tpr)
            rocAll.append(roc_auc)
            # method I: plt
            #plt.title('Receiver Operating Characteristic')
            #plt.plot(fpr, tpr, 'dodgerblue')
            
            #plt.plot([0, 1], [0, 1],'r--')
            #plt.xlim([0, 1])
            #plt.ylim([0, 1])
            #plt.ylabel('True Positive Rate')
            #plt.xlabel('False Positive Rate')

            # evaluate predictions
            scoreList = metrics.balanced_accuracy_score(y_test, predProbDf['PredictedLabel'])
            scoresXGB.append(scoreList)
            predProbDf1.append(predProbDf)
            roc_auc1.append(roc_auc)
            #CalibrationDisplay.from_predictions(y_test.values, preds, n_bins=15)
        
        #topFeatureToPlot = [X.columns.get_loc(c) for c in top20[0].tolist() if c in X]
        #shap.summary_plot(np.array(SHAP_values_per_fold)[:,topFeatureToPlot], np.array(test_values_per_fold)[:,topFeatureToPlot], feature_names = top20[0])
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        #print(len(np.mean(featureImportance, axis=0)))
        #print(len(keys))
        dataImportance = pd.DataFrame(data=np.mean(featureImportance, axis=0), index=X.columns, columns=["score"]).sort_values(by = "score", ascending=False)

        #legend_elements = [Line2D([0], [0], color='dodgerblue', lw=4, label='Mean AUC = %0.2f' % np.mean(roc_auc))]
        #plt.legend(handles=legend_elements, loc='lower right', fancybox=True, framealpha=0.5)
        #plt.savefig('roc.png', transparent=True, bbox_inches = 'tight')
        #plt.show()
        df = pd.concat(predProbDf1)
        mean1.append(np.mean(scoresXGB, axis = 0))
        std1.append(np.std(scoresXGB, axis = 0))
    print(np.mean(mean1, axis = 0),np.mean(std1, axis = 0))
    return(np.mean(mean1, axis = 0),np.mean(std1, axis = 0), df, dataImportance)



# %%
