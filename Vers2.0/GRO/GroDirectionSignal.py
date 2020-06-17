
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns

import matplotlib
import matplotlib.gridspec as gridspec

import pickle

g="A1"
A1 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])
g="A2"
A2 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


g="P1"
P1 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])
g="P2"
P2 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


g="M1"
M1 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])
g="M2"
M2 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


g="B1"
B1 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])
g="B2"
B2 = pd.read_csv("creCtr"+g+".tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])



def concat12(A1, A2, filter=False):
    A = pd.concat([A1,A2], axis=1)
    # A["dmso.-"] = (A["1dmso.-"]+1) / (A["2dmso.-"]+1)
    # A["dmso.+"] = (A["2dmso.+"]+1) / (A["1dmso.+"]+1)
    # A["dht.-"] = (A["1dht.-"]+1) / (A["2dht.-"]+1)
    # A["dht.+"] = (A["2dht.+"]+1) / (A["1dht.+"]+1)
    A["dmso.-"] = (A["1dmso.-"]+1)# / (A["A2dmso.-"]+1)
    A["dmso.+"] = (A["2dmso.+"]+1)# / (A["A1dmso.+"]+1)
    A["dht.-"] = (A["1dht.-"]+1) #/ (A["A2dht.-"]+1)
    A["dht.+"] = (A["2dht.+"]+1) #/ (A["A1dht.+"]+1)
    #
    #
    A["Directionality"] = A["dht.+"] / ( A["dht.+"] + A["dht.-"] )
    #
    if filter:
        A = A[((A["dmso.+"] > 1) & (A["dmso.-"] > 1)) | ((A["dht.-"] > 1) & (A["dht.+"] > 1))]
        # A = A[(A["Directionality"] > 0.1) & (A["Directionality"] < 0.9)]
    return A



A = concat12(A1,A2)
A["class"] = "A"
B = concat12(B1,B2)
B["class"] = "B"
P = concat12(P1,P2)
P["class"] = "P"
M = concat12(M1,M2)
M["class"] = "M"

DF = pd.concat([A,B,P,M])


classes = ("A", "B", "P", "M")

fig = plt.figure(figsize=[9,6])
gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)

fig.add_subplot(gs[0])
for cl in classes:
    p = DF[DF["class"] == cl]["Directionality"]
    sns.distplot(p, rug=False, hist=False, label=cl)
    plt.ylabel("PDF")

plt.legend()

fig.add_subplot(gs[1])
kwargs = {'cumulative': True}
for cl in classes:
    p = DF[DF["class"] == cl]["Directionality"]
    sns.distplot(p, kde_kws=kwargs, hist=False, label=cl)
    plt.ylabel("CDF")

plt.legend()

fig.savefig("directional.pdf")

pickle.dump(DF, open("Gro.DF.p", "wb"))


A = concat12(A1,A2, filter=True)
B = concat12(B1,B2, filter=True)
P = concat12(P1,P2, filter=True)
M = concat12(M1,M2, filter=True)


A["class"] = "A"
B["class"] = "B"
P["class"] = "P"
M["class"] = "M"

DF = pd.concat([A,B,P,M])

classes = ("A", "B", "P", "M")

fig = plt.figure(figsize=[9,6])
gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)

fig.add_subplot(gs[0])
for cl in classes:
    p = DF[DF["class"] == cl]["Directionality"]
    sns.distplot(p, rug=False, hist=False, label=cl)
    plt.ylabel("PDF")

plt.legend()

fig.add_subplot(gs[1])
kwargs = {'cumulative': True}
for cl in classes:
    p = DF[DF["class"] == cl]["Directionality"]
    sns.distplot(p, kde_kws=kwargs, hist=False, label=cl)
    plt.ylabel("CDF")

plt.legend()

fig.savefig("directionalFilt.pdf")
