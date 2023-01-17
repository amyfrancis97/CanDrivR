
ATAC = pd.read_csv("absSumit.bed", sep = "\t", header=None)

ATAC[1] = ATAC[1].astype('int')
ATAC[2] = ATAC[2].astype('int')

for i in [3, 4, 5, 6]:
    ATAC.iloc[:, i] = ATAC.iloc[:, i].astype('string')

# switches the regions around where pos2 > pos1, doesnt really matter about strandedness for this
x = ATAC[ATAC[1] > ATAC[2]][2]
y = ATAC[ATAC[1] > ATAC[2]][1]

indicesToChange = ATAC[ATAC[1] > ATAC[2]].index.tolist()
ATAC.loc[indicesToChange, 1] = x.tolist()
ATAC.loc[indicesToChange, 2] = y.tolist()

ATAC.to_csv("absSumitInt.bed", index=None, header = None, sep = "\t")
############### varmap #######################
varMapPrecomputed = pd.read_table("/Users/uw20204/Downloads/resultsVarMap.tsv",encoding="ISO-8859-1", on_bad_lines='skip')

varMapPrecomputed = varMapPrecomputed.drop((varMapPrecomputed[varMapPrecomputed['CHROMOSOME'].isnull()].index.tolist()))
varMapPrecomputed = varMapPrecomputed.drop((varMapPrecomputed[varMapPrecomputed['CHROMOSOME'] == 'MT'].index.tolist()))
varMapPrecomputed = varMapPrecomputed.drop((varMapPrecomputed[varMapPrecomputed['CHROMOSOME'] == 'X'].index.tolist()))
varMapPrecomputed = varMapPrecomputed.drop((varMapPrecomputed[varMapPrecomputed['CHROMOSOME'] == 'Y'].index.tolist()))

vepIDCol = "chr" + varMapPrecomputed['CHROMOSOME'].astype(int).astype('string') + "_" + varMapPrecomputed['COORDS'].astype('string') + "_" + varMapPrecomputed['USER_BASE'] + "/" + varMapPrecomputed['USER_VARIANT']
varMapPrecomputed['vepID'] = vepIDCol

# convert columns into categorical
for i in ['CAT_RES', 'DISULPHIDE', 'NTO_DNA', 'NTO_LIGAND', 'NTO_METAL', 'NTO_PROTEIN']:
    varMapPrecomputed[i] = pd.Categorical(varMapPrecomputed[i]).codes
#%%
TFBS = pd.read_table("/Users/uw20204/Downloads/TFBSencodeHG38/TFBS.bed", header = None)
TFBS['vepID'] = TFBS[0] + "_" + TFBS[1].astype('string') + "_" + TFBS[4] + "/" + TFBS[5]
TFBS = TFBS.drop([18, 17, 16, 14, 13, 12, 6, 5, 4, 3, 2, 1, 0], axis = 1)

#%%
######### green DB #########
greenDB = pd.read_csv("greenDBIntersect.txt", sep = "\t", header = None)
greenDB = greenDB.iloc[:, [0, 1, 2, 3, 4, 5, 6,9, 10, 14, 17, 21, 25]]
greenDB['vepID'] = greenDB[0] + "_" + greenDB[1].astype('string') + "_" + greenDB[3] + "/" + greenDB[4]
greenDB[greenDB["vepID"] == greenDB["vepID"].unique().tolist()[1]]
#%%
for variant in range(0, 5):
    print(greenDB["vepID"].unique().tolist()[1])
    for i in greenDB[10].unique():
        print(len(greenDB[greenDB["vepID"] == greenDB["vepID"].unique().tolist()[0]][10][greenDB[greenDB["vepID"] == greenDB["vepID"].unique().tolist()[1]][10] == i]))
#%%
for i in greenDB[10].unique():
    print(len(greenDB[greenDB["vepID"] == greenDB["vepID"].unique().tolist()[0]][10][greenDB[greenDB["vepID"] == greenDB["vepID"].unique().tolist()[1]][10] == i]))
#################