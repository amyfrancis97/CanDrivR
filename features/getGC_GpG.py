from Bio import SeqIO
import os
import pandas as pd
from strkernel.mismatch_kernel import MismatchKernel
from strkernel.mismatch_kernel import preprocess
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import classification_report # classfication summary
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import sys

if __name__ == "__main__":
    variantType = sys.argv[1]
    print(variantType)
    variants = "/bp1/mrcieu1/data/encode/public/CanDrivR/training/" + variantType + "/cosmicGnomadVariants_" + variantType + ".bed"
    outputDir = "/bp1/mrcieu1/data/encode/public/CanDrivR/training/" + variantType + "/features/"

    # reads in the human GRCh38 genome in fasta format
    record_dict = SeqIO.to_dict(SeqIO.parse("hg38_seq.fa", "fasta"))

    # reading in the variant file
    variants = pd.read_csv(variants, sep = "\t", names = ['chrom', 'pos', 'pos2', 'ref_allele', 'alt_allele', 'R','driver_status'])

    # Drops variants on the sex chromosomes
    variants = variants[(variants['chrom'] != "chrX") & (variants['chrom'] != "chrY")]
    variants = variants.reset_index(drop = True)

    # For a given variant and window size, this function queries the HG38 genome
    # Extracts the wild type sequence
    # Calculates the GC content and CpG Count
    def getGCContent(variantIndex, windowSize):
        wildType = str(record_dict[variants.loc[variantIndex, "chrom"]].seq[int(variants.loc[variantIndex, "pos"]-1-windowSize):int(variants.loc[variantIndex, "pos"]-1+windowSize)]).upper()
        GCContent = (wildType.count("G") + wildType.count("C"))/ len(wildType)
        CpGCount = (wildType.count("CG"))
        if wildType.count("G") != 0 & wildType.count("C") != 0:
            CpG_obs_exp = (CpGCount * len(wildType))/ (wildType.count("G") * wildType.count("C"))
        else:
            CpG_obs_exp = 0
        return GCContent, CpGCount, CpG_obs_exp

    # Carries out above function and appends into a dataframe
    def getGC(windowSize):
        GCDf = [getGCContent(x, windowSize) for x in range(0, len(variants))]
        GCdataframe = pd.concat([variants, pd.DataFrame(GCDf)], axis = 1)
        GCdataframe["chrom"] = GCdataframe["chrom"].str.replace("chr", "").astype(str)
        GCdataframe["pos"] = GCdataframe["pos"].astype(int)
        GCdataframe = GCdataframe.rename(columns = {0:str(windowSize*2) + "GCContent", 1:str(windowSize*2) + "CpGCount", 2:str(windowSize*2) + "CpG_obs_exp"})
        GCdataframe = GCdataframe.drop("pos1", axis = 1)
        GCdataframe.to_csv(outputDir + str(windowSize*2) + "_GC.txt", index=None, sep = "\t")

    # Applies above function to different window sizes
    for i in [50, 100, 250, 500, 1000, 5000]:
        getGC(i)
