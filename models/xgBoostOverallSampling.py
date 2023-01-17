#%%
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
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split
from sklearn.calibration import calibration_curve
import seaborn as sns

# samples positive dataset to match number of negatives
# may have to repeat this step to train on all positive data?
import warnings
warnings.filterwarnings('ignore') # setting ignore as a parameter

def testLOCO(df):
    # randomly sample dataset for balances classes
    dataframe = pd.concat([df[df.driver_status == 1].sample(len(df[df.driver_status == 0])), df[df.driver_status == 0]]).reset_index(drop=True)
    print(dataframe.groupby(['chrom', 'driver_status']).size())
    X = dataframe.drop(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"], axis=1)
    y = dataframe["driver_status"]
    groups = dataframe["grouping"]
    logo = LeaveOneGroupOut()
    predProbDf1 = []
    indices = []
    scoresXGBMean = []
    scoresXGBStd = []
    scoresXGB = []
    fprMean = []
    tprMean = []
    fprMeanFin = []
    tprMeanFin = []
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

        # Fit the model with training data and target values
        xgb_model.fit(X_train, y_train)
        y_pred = xgb_model.predict(X_test)

        predProb = xgb_model.predict_proba(X_test)
        predProbDf = pd.DataFrame(predProb)
        predProbDf.insert(0, 'PredictedLabel', y_pred)
        predProbDf.insert(0, 'ActualLabel', y_test.values)

        probs = xgb_model.predict_proba(X_test)

        # evaluate predictions
        scoreList = [metrics.accuracy_score(y_test, y_pred),
                        metrics.precision_score(y_test, y_pred)]
        scoresXGB.append(scoreList)
        predProbDf1.append(predProbDf)

    scoresXGBMean.append(np.mean(scoresXGB, axis = 0))
    scoresXGBStd.append(np.std(scoresXGB, axis = 0))
    df = pd.concat(predProbDf1)
    print(scoresXGBMean)
    return(scoresXGBMean,scoresXGBStd, df)
#%%
res = testLOCO(trainingData)
confidenceIntervals = res[2]
#%%
res
#%%
# Get ROC
preds = probs[:,1]
fpr, tpr, threshold = metrics.roc_curve(y_test, preds)

roc_auc = metrics.auc(fpr, tpr)

# method I: plt
plt.title('Receiver Operating Characteristic')
plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()
#%%
pd.DataFrame(res[0], columns = ['acc', 'precision'])
#%%
confidenceIntervals
#%%
negatives = confidenceIntervals[confidenceIntervals['PredictedLabel'] == 0][0].tolist() + confidenceIntervals[confidenceIntervals['PredictedLabel'] == 0][1].tolist()
positives = confidenceIntervals[confidenceIntervals['PredictedLabel'] == 1][0].tolist() + confidenceIntervals[confidenceIntervals['PredictedLabel'] == 1][1].tolist()
#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
bins = np.linspace(0, 1, 100)
ax.hist(negatives, bins=bins, density=True, histtype='bar', label='Negative', alpha = 0.7, color='darkblue', edgecolor='black', linewidth=1.1)
ax.hist(positives, bins=bins, density=True, histtype='bar', label='Positive', alpha = 0.7, color='lightblue', edgecolor='black', linewidth=1.1)
ax.set_xlabel('p-value')
ax.set_ylabel('density')
fig.tight_layout()
ax.legend(loc='best')
plt.show()
#%%
import pandas as pd
import numpy as np
from sklearn import datasets
 
import matplotlib.pyplot as plt
import matplotlib
%matplotlib inline
 
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['figure.dpi'] = 100
 
from IPython.core.pylabtools import figsize
figsize(8, 5)
fig, ax = plt.subplots()
bins = np.linspace(0, 1, 50)
colors=['magenta', 'turquoise']

plt.hist([negatives, positives], bins=bins, range=(0, 1), align = 'mid', density=True, histtype='bar', label=['Negative', 'Positive'], alpha = 0.7, color=colors, edgecolor='black', linewidth=1.1)
ax.set_xlabel('p-value')
ax.set_ylabel('density')
fig.tight_layout()
ax.legend(loc='best')

#

#%%
os.chdir('/Users/uw20204/Desktop/20221110/')
trainingData = pd.read_csv("featuresAll.txt", sep = "\t")
#%%


#%%
res[1]
#%%
# Run model for individual features and print output of LOCO
features = list(set(list(trainingData.columns)) - set(["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"]))
featureScoresres = []
def getScoresPerFeature(feature):
    print(feature)
    singleDf = trainingData[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping", feature]]
    return [feature, testLOCO(singleDf)[0]]
featureScores = [getScoresPerFeature(feature) for feature in features]
featureScores = pd.DataFrame(featureScores)
#%%
featureScores
#%%
dfRes = []
dfFeatures = []
for i in range(0, len(featureScores)):
    dfFeatures.append(featureScores.loc[0, i][0])
    dfRes.append(list(featureScores.loc[0, i][1]))
dfFinal = pd.DataFrame(dfRes)
dfFinal.insert(0, 'features', dfFeatures)
#%%


#%%
singleDf = trainingData[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"]]
#%%

#%%
a_list= []
resAll = []
for feature in top_features:
    a_list.append(feature)
    dfTopFeatures = df[["chrom", "pos", "driver_status", "ref_allele", "alt_allele", "reccurance", "vepID", "grouping"] + a_list]
    res = testLOCO(dfTopFeatures, 1000)
    resAll.append(res)
dfFinal = pd.DataFrame(dfRes)