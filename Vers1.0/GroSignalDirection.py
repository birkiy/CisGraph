
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns

import matplotlib
import matplotlib.gridspec as gridspec

import pickle


con1 = pd.read_csv("con.FirstPart.tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])

con2 = pd.read_csv("con.SecondPart.tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


ind1 = pd.read_csv("ind.FirstPart.tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])

ind2 = pd.read_csv("ind.SecondPart.tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


non1 = pd.read_csv("non.FirstPart.tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])

non2 = pd.read_csv("non.SecondPart.tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


tsP1 = pd.read_csv("tsP.FirstPart.tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])

tsP2 = pd.read_csv("tsP.SecondPart.tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])


tsM1 = pd.read_csv("tsM.FirstPart.tab", sep="\t", names=["1dmso.-","1dmso.+","1dht.-","1dht.+"])

tsM2 = pd.read_csv("tsM.SecondPart.tab", sep="\t", names=["2dmso.-","2dmso.+","2dht.-","2dht.+"])






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



con = concat12(con1,con2)
con["class"] = "con"
ind = concat12(ind1,ind2)
ind["class"] = "ind"
non = concat12(non1,non2)
non["class"] = "non"
tsP = concat12(tsP1,tsP2)
tsP["class"] = "tsP"
tsM = concat12(tsM1,tsM2)
tsM["class"] = "tsM"



DF = pd.concat([con, ind, non, tsP, tsM])

pickle.dump(DF, open("Gro.DF.p", "wb"))
