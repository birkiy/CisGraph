




import pandas as pd

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from statannot import add_stat_annotation

from sklearn.cluster import KMeans



metadata = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/rawData/DHS_Index_and_Vocabulary_metadata.tsv")
samples = metadata[metadata["System"] == "Genitourinary"]["DCC File ID"].values

indexARBS = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ARBS.index.txt", names=["index"])

binarydata = np.loadtxt("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/indexed.binary.txt", delimiter="\t")



systems = metadata.System.unique()[:-1]
systemDict = {}
for system in systems:
    idxSamples = metadata[metadata["System"] == system].index
    subMat = np.array([binarydata[:,i] for i in idxSamples])
    systemDict[system] = (subMat.sum(axis=0) / len(idxSamples)) * 100

systemDF = pd.DataFrame(systemDict, index=indexARBS["index"])

idxSamples = metadata[metadata["System"] == "Genitourinary"].index
subMat = np.array([binarydata[:,i] for i in idxSamples])

binDF = pd.DataFrame(data=np.transpose(subMat), columns=samples, index=indexARBS["index"])




header = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "index.ARBS"]
con = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/con.cre.bed", names=header)
con["nodeClass"] = "con"
ind = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ind.cre.bed", names=header)
ind["nodeClass"] = "ind"
non = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/non.cre.bed", names=header)
non["nodeClass"] = "non"

cre = pd.concat([con, ind, non]).reset_index(drop=True)
cre["nlog"] = np.log(cre["numsamples"])

systemDict = systemDF.to_dict(orient="index")
binDict = binDF.to_dict(orient="index")
creDict = cre.to_dict(orient="index")



concatDict = {}
for i in range(cre.shape[0]):
    idx = cre.loc[i, "index.ARBS"]
    tmp = dict(**creDict[i], **binDict[idx], **systemDict[idx])
    concatDict[i] = tmp

creBin = pd.DataFrame.from_dict(concatDict, orient="index")

creBin = creBin.sort_values(["nodeClass", "Genitourinary"], ascending=False)




Nc = range(1, 20)
kmeans = [KMeans(n_clusters=i) for i in Nc]
score = [kmeans[i].fit(creBin[systems]).score(creBin[systems]) for i in range(len(kmeans))]


kmeans = KMeans(n_clusters=4, random_state=0).fit(creBin[systems])

creBin['cluster'] = kmeans.labels_

creBin.groupby(["nodeClass", "cluster"]).size()


############################
# PCA

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
features = systems

x = creBin.loc[:, features].values

y = creBin.loc[:, ["nodeClass"]].values

x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])


finalDf = pd.concat([principalDf, creBin[['nodeClass']]], axis = 1)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = ['con', 'ind', 'non']
# colors = ["#63b7af", "#abf0e9", "#d4f3ef"]
colors = "rgb"
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['nodeClass'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)


ax.legend(targets)
ax.grid()


fig.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/tissues.DHS.pdf")


##################################3

msk = np.random.rand(len(creBin)) < 0.8
train = creBin[msk]
test = creBin[~msk]


y = np.array([1 if i == "con" else 2 if i == "ind" else 4 if i == "non" else None for i in train["nodeClass"]])

m = 3
y = np.transpose((((y[:,None] & (1 << np.arange(m)))) > 0).astype(int))
x = np.transpose(train[systems].values) / 100




def acc(w, x, y):
    return np.mean([i == k for i, k in zip(np.argmax(np.dot(w,x), axis=0),np.argmax(y, axis=0))])



def softmax(w,x,y,lr=0.01):
    probs = np.exp(np.dot(w,x))
    probs /= np.sum(probs)
    error = probs - y
    w -= lr * np.tensordot(error, x, axes=0)

import random

def perceptron(w,x,y):
    gu = np.argmax(np.dot(w,x))
    cl = np.argmax(y)
    if gu != cl:
        w[cl,:] += x
        w[gu,:] -= x



def train(algo, x,y,T=2**20):
    nTrn, xDim, yDim = x.shape[1], x.shape[0], y.shape[0]
    w = np.random.rand(yDim, xDim)
    nextPrint = 1
    for t in range(T):
        i = random.randint(0,nTrn-1)
        algo(w, x[:,i], y[:,i])
        if t == nextPrint:
            print(f"{i}th sample with {acc(w,x,y)} accuracy.")
            nextPrint = min((2*t, T))
    return w

train(softmax, x, y)


###########################

# Heatmap

lut = dict(zip(creBin["nodeClass"].unique(), ["#63b7af", "#abf0e9", "#d4f3ef"]))
row_colors = creBin["nodeClass"].map(lut)

g = sns.clustermap(creBin[systems], row_colors=row_colors, cmap="YlOrBr", row_cluster=False)

g.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/heatmapTissues.pdf")


######################################3

prostateTissues = [
    'ENCFF503PAE', 'ENCFF782HVG',
    'ENCFF939VCZ', 'ENCFF329UDK', 'ENCFF015ICH', 'ENCFF658UUI',
    'ENCFF335MSP', 'ENCFF537ICV', 'ENCFF081ETZ', 'ENCFF875DXE',
    'ENCFF898KEE', 'ENCFF309DYE', 'ENCFF176OMF', 'ENCFF625QKW'
    ]
labels = ['Hela', 'HeLaS3', 'PrEC', 'PrEC', 'LNCap', 'LNCap', 'fTestes',
       'fTestes', 'Ovary', 'fOvary', 'HeLaS3', 'HeLaS3', 'PC3', 'PC3']

creBin = creBin.rename(columns=dict(zip(prostateTissues, labels)))


creBin[prostateTissues]


lut = dict(zip(creBin["nodeClass"].unique(), ["#63b7af", "#abf0e9", "#d4f3ef"]))
row_colors = creBin["nodeClass"].map(lut)


g = sns.clustermap(creBin[labels], row_colors=row_colors, cmap="YlOrBr", row_cluster=False)

g.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ProstateTissues.pdf")



creDF = pd.concat([cre, binDF])


binarydata = np.loadtxt("/groups/lackgrp/ll_members/berkay/DHS/rawData/dat_bin_FDR01_hg19.txt", delimiter="\t")


file = "/groups/lackgrp/ll_members/berkay/DHS/rawData/dat_bin_FDR01_hg19.txt"




import csv
bed = {}
with open(file) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for i, row in enumerate(reader):
        # row3 = row[3].split(".")[0]
        row3 = [row[i] for i in genTissues]
        bed[i] = row3
