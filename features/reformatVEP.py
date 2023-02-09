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
testing = [getAA(variant) for variant in range(0, len(df))]
df = df.fillna(0)
df = df.drop([2, 6, 7], axis =1 )
df = df.rename(columns = {0:'chrom', 1: 'pos', 3: 'ref_allele', 4: 'alt_allele', 5:'driver_stat'})

# Save table
df.to_csv(featureOutputDir + "vepAA.bed", index=None, sep = "\t")

