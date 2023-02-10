# Import packages
import pandas as pd
import os
import re
from textwrap import wrap
import numpy as np
from config import *

################### Get Consequences from VEP ###################
# Read in the VEP ouput file
df = pd.read_csv(vepConsequenceOutput, sep = "\t", header =None, skiprows=5)

# Create one hot autoencoding of consequences from VEP
def getConsequences(variant):
    consList = list(set(re.split(r'=|,|&', df[7][variant])[1:]))
    df.loc[variant, consList] = 1
testing = [getConsequences(variant) for variant in range(0, len(df))]
df = df.fillna(0)
df = df.drop([2, 6, 7], axis =1 )
df = df.rename(columns = {0:'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5:'driver_stat'})

# Save table
df.to_csv(featureOutputDir + "vepConsequences.bed", index=None, sep = "\t")

################# Get Amino Acids from VEP ##################
# Read in the VEP ouput file
df = pd.read_csv(vepAAOutput, sep = "\t", header =None, skiprows=5)

def getAA(variant):
    WT = re.split(r'=|,|/', df[7][variant])[1]
    if WT == "": # if there are no amino acids available then the variant is located in a non-protein-coding region
        aminoAcids="none"
        WT = "-"
        mutant = "-"
    else:
        if len(re.split(r'=|,|/', df[7][variant])) == 2:
            mutant = WT
        else:
            mutant = re.split(r'=|,|/', df[7][variant])[2]
            if mutant == "": # If there is no mutant then change is synonymous
                mutant = WT
            else:
                mutant = mutant
            # Reformat to one-hot-encoding
            df.loc[variant, "WT_" + WT] = 1
            df.loc[variant, "mutant_" + mutant] = 1
    return[WT, mutant]
testing = [getAA(variant) for variant in range(0, len(df))]
df = df.fillna(0)
df = df.drop([2, 6, 7], axis =1 )
df = df.rename(columns = {0:'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5:'driver_stat'})

# Get amino acids in a non-one-hot-autoencoded format
AA = pd.DataFrame(testing)
AA = AA.rename(columns = {0:'WT_AA', 1:'mutant_AA'})
AA = pd.concat([df[['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat']], AA], axis = 1)

# Save tables
df.to_csv(featureOutputDir + "vepAA_OHA.bed", index=None, sep = "\t")
AA.to_csv(featureOutputDir + "vepAA.bed", index=None, sep = "\t")


############### Get distances to closest feature from VEP ################
# Read in the VEP ouput file
df = pd.read_csv(vepDistanceOutput , sep = "\t", header =None, skiprows=5)

# Create one hot autoencoding of consequences from VEP
def getDistance(variant):
    distanceList = list(set(re.split(r'=|,|&', df[7][variant])[1:]))
    if distanceList != [""]:
        distanceList = list(filter(None, distanceList))
        distanceList = [int(distance) for distance in distanceList]
        meanDistance = round(np.mean(distanceList), 1)
        minDistance = round(np.min(distanceList), 1)
        df.loc[variant, "meanDistance"] = meanDistance
        df.loc[variant, "minDistance"] = minDistance
    else:
        meanDistance = []
        minDistance = []
    return [meanDistance, minDistance]
testing = [getDistance(variant) for variant in range(0, len(df))]

df = df.fillna(0)
df = df.drop([2, 6, 7], axis =1 )
df = df.rename(columns = {0:'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5:'driver_stat'})

# Save table
df.to_csv(featureOutputDir + "vepDistance.bed", index=None, sep = "\t")
