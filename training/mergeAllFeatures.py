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
dnaShape = dnaShape.drop(['end', 'driver_stat'], axis = 1)
dnaShape = dnaShape.rename(columns = {'start': 'pos', 'ref' : 'ref_allele', 'alt': 'alt_allele'})
dnaShape['vepID'] = dnaShape['chrom'] + "_" + dnaShape['pos'].astype('string') + "_" + dnaShape['ref_allele'] + "/" + dnaShape['alt_allele'] 

#%%
LS_annotation = pd.read_csv("LS_annotation.txt", sep = "\t")
#%%
AS_annotation = pd.read_csv("AS_annotation.txt", sep = "\t")
AS_annotation = AS_annotation.drop("pos1")
AS_annotation = AS_annotation.rename(columns = {'pos2': 'pos'})
#%%
dinucleotideProperties = pd.read_csv("dinucleotideProperties.txt", sep = ",")
dinucleotideProperties = dinucleotideProperties.rename(columns = {'start': 'pos'}).drop(['end', 'WTtrinuc', 'mutTrinuc'], axis =1)
#%%
dnaShape
#%%
# Merge all of the current feature datasets together
df = cons.merge(vep)
df = df.merge(vepAA)
df = df.merge(dfKernel)
df = df.merge(dfGC)
df = df.merge(AAproperties)
df = df.merge(AASubstMatrix)
df = df.merge(ATAC2)
df = df.merge(TSSdistance)
df = df.merge(dfUnique)
#%%
df = df.merge(dnaShape, on = 'vepID')
#%%df = df.merge(LS_annotation)
#%%
#df = df.merge(AS_annotation)
df = df.merge(dinucleotideProperties)

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
