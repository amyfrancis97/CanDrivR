# Reformat VEP output
import pandas as pd
import ast
import re
import os
os.chdir('/Users/uw20204/Desktop/20221110')

# Read in Cosmic & Gnomad VEP tables
cosmicVEP1 = pd.read_csv("cosmicVEP.txt", sep = "\t")
gnomadVEP1 = pd.read_csv("gnomadVEP.txt", sep = "\t")

# Get amino acid changes from VEP
def getVEPaa(dataset):
    AA = pd.concat([dataset.drop_duplicates('#Uploaded_variation', keep='first')['#Uploaded_variation'], dataset.drop_duplicates('#Uploaded_variation', keep='first')['Amino_acids'].str.split("/", expand = True)], axis = 1).reset_index(drop = True)
    AA = AA.rename(columns = {'#Uploaded_variation': "vepID", 0:"AA1", 1:"AA2"})
    AA["AA2"][AA["AA2"].isnull()] = list(AA[AA["AA2"].isnull()]["AA1"])
    AA[""] = pd.factorize(AA["AA2"], sort=True)[0] + 1 
    AA["AA1"] = pd.factorize(AA["AA1"], sort=True)[0] + 1
    return AA

vepAACosmic = getVEPaa(cosmicVEP1)
vepAAGnomad= getVEPaa(gnomadVEP1)

# assigning driver status label before concatenating
vepAACosmic['driver_status'] = 1
vepAAGnomad['driver_status'] = 0
vepAA = pd.concat([vepAACosmic, vepAAGnomad])

# Get consequences from VEP
def getVEP(dataset):
    groups = dataset.groupby('#Uploaded_variation')['Consequence'].apply(set).apply(list)
    df1 = groups.reset_index(name 
                            = 'listvalues')

    consequences = list(set(list(dataset['Consequence'].str.split(",", expand = True).stack())))
    newdf = pd.DataFrame(columns = ['variant'] + consequences)
    newdf['variant'] = df1['#Uploaded_variation']
    
    # Converting VEP format to binary
    def fillConsDf(i):

        newdf.iloc[i, 1:] = 0 # Automatically fill all consequences with '0'
        listtoflatten = [df1.iloc[i, 1][x].split(",") for x in range(0, len(df1.iloc[i, 1]))]
        flattened = [val for sublist in listtoflatten for val in sublist]
        indexCols = list(set(flattened))
        newdf.loc[i, indexCols] = 1 # Replace consequence columns present with a '1'

    output = [fillConsDf(i) for i in range(0,len(newdf))]
    return(newdf)

gnomadVEP = getVEP(gnomadVEP1)
cosmicVEP = getVEP(cosmicVEP1)

# assigning driver status label before concatenating
cosmicVEP['driver_status'] = 1
gnomadVEP['driver_status'] = 0

vep = pd.concat([gnomadVEP, cosmicVEP])
vep = vep.rename(columns = {'variant':'vepID'})

# Write to csv
vep.to_csv("featuresTraining/vep.txt", index=None, header = None, sep = "\t")
vepAA.to_csv("featuresTraining/vepAA.txt", index=None, header = None, sep = "\t")
