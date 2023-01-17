#%%
import pandas as pd
import ast
import re
import os
# xgboost packages
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

# samples positive dataset to match number of negatives
# may have to repeat this step to train on all positive data?
import warnings
warnings.filterwarnings('ignore') # setting ignore as a parameter

def testLOCO(df, TrainSampleSize):
    scoresXGB3 = []
    for i in range(0,5):
        sampled = pd.concat([df[df.driver_status == 1].sample(len(df[df.driver_status == 0])), df[df.driver_status == 0]]).reset_index(drop=True)
        dataframe = sampled

        #print(dataframe.groupby('chrom').count()['pos'])

        #print(sampled.groupby('chrom').count())

        X = dataframe.drop(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"], axis=1)
        y = dataframe["driver_status"]
        groups = dataframe["grouping"]
        logo = LeaveOneGroupOut()

        indices = []
        scoresXGB2 = []
        for train_index, test_index in logo.split(X, y, groups):
            indices.append([list(train_index), (list(test_index))])
            X_train1, X_test = X.iloc[train_index, :].reset_index(drop = True), X.iloc[test_index, :].reset_index(drop = True)
            y_train1, y_test = y[train_index].reset_index(drop = True), y[test_index].reset_index(drop = True)
            

            scoresXGB = []
            for i in range(0, 5):
                # randomly sample training dataset for balances classes
                neg = list(y_train1[y_train1 == 0].sample(TrainSampleSize).index)
                pos = list(y_train1[y_train1 == 1].sample(TrainSampleSize).index)
                balancedIndices = neg + pos
                y_train = y_train1[balancedIndices]
                X_train = X_train1.iloc[balancedIndices, :]
                from numpy import asarray
                from sklearn.preprocessing import StandardScaler
                # define standard scaler
                scaler = StandardScaler()
                # transform data
                X_train = scaler.fit_transform(X_train)
                X_test = scaler.fit_transform(X_test)

                xgb_model =  XGBClassifier(
                learning_rate =0.001, #put this back to 0.001 in the end
                n_estimators=1000,
                max_depth=3,
                min_child_weight=4,
                gamma=1,
                subsample=0.65,
                colsample_bytree=0.6,
                objective= 'binary:logistic',
                nthread=4,
                scale_pos_weight=1,
                seed=27)
                # Fit the model with training data and target values
                xgb_model.fit(X_train, y_train)
                y_pred = xgb_model.predict(X_test)

                # plot feature importance
                #xgb_model.get_booster().feature_names = list(X.columns)
                #plot_importance(xgb_model.get_booster())
                #pyplot.show()

                predictions = [round(value) for value in y_pred]

                # evaluate predictions
                scoreList = [metrics.balanced_accuracy_score(y_test, predictions), 
                                accuracy_score(y_test, predictions),
                                metrics.precision_score(y_test, predictions),
                                metrics.recall_score(y_test, predictions),
                                metrics.roc_auc_score(y_test, predictions),
                                metrics.f1_score(y_test, predictions)]
                scoresXGB.append(scoreList)
            scoresXGB2.append(np.mean(scoresXGB, axis = 0))
        #print(np.mean(scoresXGB2, axis = 0))
        scoresXGB3.append(np.mean(scoresXGB2, axis = 0))
    print(np.mean(scoresXGB3, axis = 0))
    return(np.mean(scoresXGB3, axis = 0))
#%%

os.chdir('/Users/uw20204/Desktop/20221110/')
#%%
gnomadVEP1 = pd.read_csv("gnomadVEP.txt", sep = "\t")
cosmicVEP1 = pd.read_csv("cosmicVEP.txt", sep = "\t")

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

os.chdir('/Users/uw20204/Desktop/20221110/conservation')
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

os.chdir('/Users/uw20204/Desktop/20221110')

# Merging all conservation datasets
cons = phastCons20.merge(phyloP20, on = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', 'driver_status', 'vepID'])
for consFeature in [phyloP100, phastCons100, phastCons4, phyloP4, phastCons7, phyloP7, phastCons30, phyloP30, phyloP17, phastCons17]:
    cons = cons.merge(consFeature, on = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'reccurance', 'driver_status', 'vepID'])
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

AAproperties = pd.read_csv("VEP_web_aminoAcid_properties.txt", sep = ",")
AAproperties = AAproperties.drop(['AA1', 'AA2'], axis = 1)

ATAC = pd.read_csv("intersectAbsSumit.bed", sep = "\t", header=None)
ATAC = ATAC.drop([2, 6, 7, 8], axis = 1)
ATAC = ATAC.rename(columns = {0: 'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5: 'driver_status', 9:'mean_pileup', 10:'median_pileup', 11:'mean_fold_enrichment', 12: 'median_fold_enrichment', 13: 'distance_feature_to_variant'})
ATAC['vepID'] = ATAC['chrom'] + "_" + ATAC['pos'].astype('string') + "_" + ATAC['ref_allele'] + "/" + ATAC['alt_allele'] 

# takes the mean of any duplicates
ATAC2 = ATAC.groupby('vepID').mean().reset_index()
ATAC2['chrom'] = ATAC2['vepID'].str.split("_", expand = True)[0]
ATAC2['pos'] = ATAC2['vepID'].str.split("_", expand = True)[1].astype('int')
ATAC2['ref_allele'] = ATAC2['vepID'].str.split("_", expand = True)[2].str.split("/", expand = True)[0]
ATAC2['alt_allele'] = ATAC2['vepID'].str.split("_", expand = True)[2].str.split("/", expand = True)[1]

TSSdistance = pd.read_csv("intersectTSS.bed", sep = "\t", header = None)
TSSdistance = TSSdistance.drop([2, 6, 7, 8, 9, 10, 11, 12, 13, 14], axis =1)

TSSdistance = TSSdistance.rename(columns = {0: 'chrom', 1: 'pos', 3:'ref_allele', 4: 'alt_allele', 5: 'driver_status', 15: 'distanceTSS'})

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

vepAA = pd.read_csv("featuresTraining/vepAA.txt", sep = "\t")
vep = pd.read_csv("featuresTraining/vep.txt", sep = "\t")
#cons = pd.read_csv("featuresTraining/cons.txt", sep = "\t")

AASubstMatrix = pd.read_csv("AASubstMatrices.txt", sep = ",")
AASubstMatrix = AASubstMatrix.drop(['AA1', 'AA2'], axis =1)

#%%
dnaShape = pd.read_csv("dnaShape.txt", sep = ",")
dnaShape = dnaShape.dropna(axis=1, how='all')
dnaShape = dnaShape.drop('end', axis = 1)
dnaShape = dnaShape.rename(columns = {'start': 'pos', 'ref' : 'ref_allele', 'alt': 'alt_allele', 'driver_stat' : 'driver_status'})
dnaShape['chrom'] = dnaShape['chrom'].astype('string')
#%%
LS_annotation = pd.read_csv("LS_annotation.txt", sep = "\t")
LS_annotation
#%%
BS_annotation = pd.read_csv("BS_annotation.txt", sep = "\t")
BS_annotation = BS_annotation.rename(columns = {"distancesToLinearSpliceJunctionsSum": "distancesToBSJunctionsSum"})
#%%

AS_annotation = pd.read_csv("AS_annotation.txt", sep = "\t")
#%%
AS_annotation
#%%

dinucleotideProperties = pd.read_csv("dinucleotideProperties.txt", sep = ",")
dinucleotideProperties = dinucleotideProperties.rename(columns = {'start': 'pos'}).drop(['end', 'WTtrinuc', 'mutTrinuc'], axis =1)
new_names = [(i,i+'dinucleotideProperties') for i in dinucleotideProperties.iloc[:, 5:].columns.values]
dinucleotideProperties.rename(columns = dict(new_names), inplace=True)
#%%
df
#%%
# Merge all of the current feature datasets together
df = cons.merge(vep)
df = df.merge(vepAA)
df = df.merge(dfKernel)
df = df.merge(dfGC)
df = df.merge(AAproperties)
df = df.merge(AASubstMatrix)
df = df.merge(ATAC2)
#%%
df = df.merge(TSSdistance)
df = df.merge(dfUnique)
df = df.merge(dnaShape)
#%%
df = df.merge(LS_annotation)
df = df.merge(AS_annotation)
df = df.merge(dinucleotideProperties)
df
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
#%%

#%%
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
testLOCO(df, 1500)
#%%
df
############################################
#%%
dfTop = df.drop(featureScoresFinal[featureScoresFinal.sort_values(1, ascending=False)[1] < 0.51][0].to_list(), axis = 1)
#%%
LS_annotation.columns.tolist()[6:]
#%%
singleDf = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"]+ LS_annotation.columns.tolist()[6:]]
singleDf
#%%
#################
# Run model for individual features and print output of LOCO
#%%
features = list(set(list(df.columns)) - set(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"]))
featureScoresres = []
def getScoresPerFeature(feature):
    print(feature)
    singleDf = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping", feature]]#
    featureScoresres.append([feature, testLOCO(singleDf, 400)])
    return featureScoresres

featureScores = [getScoresPerFeature(feature) for feature in features]
featureScores = pd.DataFrame(featureScores)

dfRes = []
dfFeatures = []
for i in range(0, len(featureScores)):
    dfFeatures.append(featureScores.loc[0, i][0])
    dfRes.append(list(featureScores.loc[0, i][1]))
dfFinal = pd.DataFrame(dfRes)
dfFinal.insert(0, 'features', dfFeatures)
dfFinal.sort_values(0, ascending=False)[dfFinal.sort_values(0, ascending=False)[0] > 0.5]
#%%
dfFinal.sort_values(0, ascending=False)[dfFinal.sort_values(0, ascending=False)[0] > 0.51]
#%%
# select feature groups with accuracy >0.52

#%%
top_features = dfFinal.sort_values(0, ascending=False)[dfFinal.sort_values(0, ascending=False)[0] > 0.5]['features'].to_list()
dfTopFeatures = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"] + top_features + LS_annotation.columns.tolist()[6:]]
#%%
dfTopFeatures = dfTopFeatures.fillna(0)
#%%
len(df.columns)
#%%
#top_features = featureScoresFinal.sort_values(1, ascending=False)[featureScoresFinal.sort_values(1, ascending=False)[1] > 0.5][0].to_list()

#%%
a_list= []
resAll = []
for feature in top_features:
    a_list.append(feature)
    print(a_list)
    dfTopFeatures = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"] + a_list]
    res = testLOCO(dfTopFeatures, 1000)
    resAll.append(res)
#%%
pd.DataFrame(resAll)[0].to_list()


#%%
# Save features to CSV
vepAA.to_csv("featuresTraining/vepAA.txt", sep = "\t", index=False)
vep.to_csv("featuresTraining/vep.txt", sep = "\t", index=False)
cons.to_csv("featuresTraining/cons.txt", sep = "\t", index=False)
dfGC.to_csv("featuresTraining/dfGC.txt", sep = "\t", index=False)
#%%
