import pandas as pd
import ast
import re
import os
import sys
from config import *

if __name__ == "__main__":
    alleleFreq = float(sys.argv[1])
    type = str(sys.argv[2])
    print(alleleFreq)
    print(type)
    os.chdir(str(sys.argv[3]))
    try:
        os.remove(cosmicGnomadDir + str(type) + str(alleleFreq) + "_gnomad_AF_filtered.out")
    except OSError:
        pass
    gnomadNegatives = []
    with pd.read_table(str(type) + str(alleleFreq) + "_gnomadAlleleFreqIncl.txt", on_bad_lines='skip', header= None, chunksize=100000, engine='python') as reader:
        for chunk in reader:
            chunk = chunk.reset_index(drop=True)
            negative = []

            def getAlleleFreq(i):
                if 'controls_AF_raw=' in chunk[7][i]:
                    if float(chunk[7][i].split("controls_AF_raw=")[1].split(";")[0]) > alleleFreq:
                        negative.append(i)
                    elif 'non_cancer_AF=' in chunk[7][i]:
                    # pull out the allele frequency in the non-cancer population
                    # if the allele frequency is >3% them pull out the index
                        if float(chunk[7][i].split("non_cancer_AF=")[1].split(";")[0]) > alleleFreq:
                            negative.append(i)
                elif 'non_cancer_AF=' in chunk[7][i]:
                # pull out the allele frequency in the non-cancer population
                # if the allele frequency is >3% them pull out the index
                    if float(chunk[7][i].split("non_cancer_AF=")[1].split(";")[0]) > alleleFreq:
                        negative.append(i)
            [getAlleleFreq(i) for i in range(0, len(chunk))]
            output_path=str(type) + str(alleleFreq) + "_gnomad_AF_filtered.out"
            chunk.loc[negative][range(0, 6)].to_csv(output_path, mode='a', sep = "\t", index = None, header=None)
