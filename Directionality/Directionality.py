

import pandas as pd
import numpy as np

from scipy import stats
from statannot import add_stat_annotation

import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', 20)

Mat = pd.read_csv("concatCRE.tab", sep="\t")

pairs = (("coP.FirstPart.bed", "coP.SecondPart.bed"),
         ("coM.FirstPart.bed", "coM.SecondPart.bed"),
         ("inP.FirstPart.bed", "inP.SecondPart.bed"),
         ("inM.FirstPart.bed", "inM.SecondPart.bed"),
         ("noP.FirstPart.bed", "noP.SecondPart.bed"),
         ("noM.FirstPart.bed", "noM.SecondPart.bed"),
         ("tsP.FirstPart.bed", "tsP.SecondPart.bed"),
         ("tsM.FirstPart.bed", "tsM.SecondPart.bed"),
         ("con.FirstPart.bed", "con.SecondPart.bed"),
         ("ind.FirstPart.bed", "ind.SecondPart.bed"),
         ("non.FirstPart.bed", "non.SecondPart.bed"),
         ("tss.FirstPart.bed", "tss.SecondPart.bed"))

CrePair = pd.DataFrame()
for pair in pairs:
    pair1 = Mat[Mat["deepTools_group"] == pair[0]]
    pair2 = Mat[Mat["deepTools_group"] == pair[1]]
    pair1 = pair1[["#chrom", "start", "end", "name", "deepTools_group", "groseq_dht.plus", "groseq_dht.minus", "groseq_dmso.plus", "groseq_dmso.minus"]]
    pair2 = pair2[["#chrom", "start", "end", "name", "deepTools_group", "groseq_dht.plus", "groseq_dht.minus", "groseq_dmso.plus", "groseq_dmso.minus"]]
    pair1 = pair1.reset_index(drop=True)
    pair2 = pair2.reset_index(drop=True)
    oK = sum(pair1["end"] == pair2["start"])
    if len(pair1) == oK:
        print("Good 2 Go")
    else:
        print("No Good Bro!!")
    dhtP = (pair2["groseq_dht.plus"] + 1)
    dhtM = (pair1["groseq_dht.minus"] + 1)
    dmsP = (pair2["groseq_dmso.plus"] + 1)
    dmsM = (pair1["groseq_dmso.minus"] +1)
    crePair = pair1[["#chrom", "start" , "end"]]
    crePair["end"] += 350
    crePair["name"] = pair1["name"]
    crePair["nodeClass"] = pair[0].split(".")[0]
    crePair["P.DHT"] = dhtP
    crePair["M.DHT"] = dhtM
    crePair["P.DMSO"] = dmsP
    crePair["M.DMSO"] = dmsM
    crePair["D.DHT"] = dhtP / (dhtP + dhtM)
    crePair["D.DMSO"] = dmsP / (dmsP + dmsM)
    CrePair = pd.concat([CrePair, crePair])


CrePair["P.DHT.Log"] = np.log(CrePair["P.DHT"])
CrePair["M.DHT.Log"] = np.log(CrePair["M.DHT"])
CrePair["Lfc.DHT"] = np.log(CrePair["P.DHT"] / CrePair["M.DHT"])
CrePair = CrePair.sort_values(by=["nodeClass"], ascending=False)




CrePair.to_csv("Directional.BPM.bed", sep="\t", index=False)


CrePair = pd.read_csv("Directional.BPM.bed", sep="\t")
CreUni = CrePair[(CrePair["nodeClass"] == "tss") | (CrePair["nodeClass"] == "con" )| (CrePair["nodeClass"] == "ind") | (CrePair["nodeClass"] == "non")]

CreUni["Activity.DHT"] = CreUni["P.DHT"] + CreUni["M.DHT"]
CreUni["Activity.DMSO"] = CreUni["P.DMSO"] + CreUni["M.DMSO"]
CreUni["quantile"] = "biDirectional"

CreUni.loc[((CreUni["D.DHT"] < 0.45) & (CreUni["D.DHT"] > 0.1)), "quantile"] = "intMinus"

CreUni.loc[((CreUni["D.DHT"] > 0.55) & (CreUni["D.DHT"] < 0.9)), "quantile"] = "intPlus"

CreUni.loc[((CreUni["D.DHT"] > 0.9) & (CreUni["M.DHT"] == 1)), "quantile"] = "uniPlusOnly"

CreUni.loc[((CreUni["D.DHT"] < 0.1) & (CreUni["P.DHT"] == 1)), "quantile"] = "uniMinusOnly"

CreUni.loc[((CreUni["D.DHT"] > 0.9) & (CreUni["M.DHT"] != 1)), "quantile"] = "uniPlus"

CreUni.loc[((CreUni["D.DHT"] < 0.1) & (CreUni["P.DHT"] != 1)), "quantile"] = "uniMinus"

# CreUni = CreUni[((CreUni["nodeClass"] == "tss") & ((CreUni["P.DHT"] != 1) | (CreUni["M.DHT"] != 1))) | (CreUni["nodeClass"] == "con") | (CreUni["nodeClass"] == "non") | (CreUni["nodeClass"] == "ind")]
# CreUni = CreUni[~((CreUni["P.DHT"] == 1) & (CreUni["M.DHT"] == 1))]
# (CreUni["nodeClass"] == "tss") &

fig = plt.figure(figsize=[6,9])
g = sns.boxplot(x="nodeClass", y="Activity.DHT", data=CreUni, palette=colorPalette)
g.set_yscale("log")
boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]



add_stat_annotation(g, x="nodeClass", y="Activity.DHT", data=CreUni,  test="Mann-Whitney", box_pairs=boxPairs)


fig.savefig("Activity.pdf")




fig = plt.figure(figsize=[6,9])
g = sns.boxplot(x="nodeClass", y="Activity.DMSO", data=CreUni, palette=colorPalette)
g.set_yscale("log")
boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]



add_stat_annotation(g, x="nodeClass", y="Activity.DMSO", data=CreUni,  test="Mann-Whitney", box_pairs=boxPairs)


fig.savefig("Activity.DMSO.pdf")





fig = plt.figure(figsize=[20,9])

gs = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[4, 4])
#### DHT

# colorPalette = ["#ee8572", "#d4f3ef", "#abf0e9", "#63b7af"]
colorPalette = (("tss", "#ee8572"),
                ("non", "#d4f3ef"),
                ("ind", "#abf0e9"),
                ("con", "#63b7af"))

fig.add_subplot(gs[0])
for nodeClass, color in colorPalette:
    Z = CreUni[CreUni["nodeClass"] == nodeClass]["D.DHT"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)



plt.legend()

plt.ylabel("CDF")
plt.xlabel("D")
plt.title("DHT")

##### DMSO


fig.add_subplot(gs[1])

for nodeClass, color in colorPalette:
    Z = CreUni[CreUni["nodeClass"] == nodeClass]["D.DMSO"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend()

plt.ylabel("CDF")
plt.xlabel("D")
plt.title("DMSO")

fig.savefig("CDF.Directionality.pdf")



fig = plt.figure(figsize=[20,20])

gs = gridspec.GridSpec(ncols=2, nrows=2)

# Minus DHT
fig.add_subplot(gs[0,1])

for nodeClass, color in colorPalette:
    Z = 1/CreUni[CreUni["nodeClass"] == nodeClass]["M.DHT"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend()
plt.ylabel("CDF")
plt.xlabel("1/ReadDensity")
plt.title("Minus Sense - DHT")

### Plus DHT
fig.add_subplot(gs[0,0])

for nodeClass, color in colorPalette:
    Z = 1/CreUni[CreUni["nodeClass"] == nodeClass]["P.DHT"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.title("Plus Sense - DHT")
plt.legend()
plt.ylabel("CDF")
plt.xlabel("1/ReadDensity")


# Minus DMSO
fig.add_subplot(gs[1,1])
for nodeClass, color in colorPalette:
    Z = 1/CreUni[CreUni["nodeClass"] == nodeClass]["M.DMSO"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend()
plt.ylabel("CDF")
plt.xlabel("1/ReadDensity")
plt.title("Minus Sense - DMSO")



### Plus DMSO

fig.add_subplot(gs[1,0])
for nodeClass, color in colorPalette:
    Z = 1/CreUni[CreUni["nodeClass"] == nodeClass]["P.DMSO"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)


plt.title("Plus Sense - DMSO")
plt.legend()
plt.ylabel("CDF")
plt.xlabel("1/ReadDensity")

fig.savefig("OI.pdf")


# from numpy import random
#
# Z = random.exponential(scale=2, size=1000)
#
# N = len(Z)
# X2 = np.sort(Z)
# F2 = np.array(range(N))/float(N)
# plt.plot(X2, F2, label="ind")
#
# plt.legend()
#
# plt.ylabel("CDF")
# plt.xlabel("D")
#
# plt.show()

















#
#
# pd.cut(CreUni["D.DHT"], q=[0,0.1, 0.25, 0.75, 0.9, 1])
#
#
#
# g = sns.catplot(y="D.DHT", x="quantile",hue="nodeClass", data=CreUni, kind="violin")
# boxPairs = [(("intHigh", "tss"), ("intHigh", "con")),
#             (("intHigh", "tss"), ("intHigh", "ind")),
#             (("intHigh", "tss"), ("intHigh", "non")),
#             (("intHigh", "con"), ("intHigh", "ind")),
#             (("intHigh", "con"), ("intHigh", "non")),
#             (("intHigh", "ind"), ("intHigh", "non"))]
#
#
# add_stat_annotation(g.fig.get_axes()[0], y="D.DHT", x="quantile",hue="nodeClass", data=CreUni,  test="Mann-Whitney", box_pairs=boxPairs)
#
#


colorPalette = ["#ee8572", "#d4f3ef", "#abf0e9", "#63b7af"]

fig = plt.figure(figsize=[12,9])
g = sns.boxplot(x="quantile", y="Activity", hue="nodeClass", data=CreUni, palette=colorPalette)
g.set_yscale("log")
boxPairs = [(("intPlus", "tss"), ("intPlus", "con")),
            (("intPlus", "tss"), ("intPlus", "ind")),
            (("intPlus", "tss"), ("intPlus", "non")),
            (("intPlus", "con"), ("intPlus", "ind")),
            (("intPlus", "con"), ("intPlus", "non")),
            (("intPlus", "ind"), ("intPlus", "non")),
            (("biDirectional", "tss"), ("biDirectional", "con")),
            (("biDirectional", "tss"), ("biDirectional", "ind")),
            (("biDirectional", "tss"), ("biDirectional", "non")),
            (("biDirectional", "con"), ("biDirectional", "ind")),
            (("biDirectional", "con"), ("biDirectional", "non")),
            (("biDirectional", "ind"), ("biDirectional", "non")),
            (("intMinus", "tss"), ("intMinus", "con")),
            (("intMinus", "tss"), ("intMinus", "ind")),
            (("intMinus", "tss"), ("intMinus", "non")),
            (("intMinus", "con"), ("intMinus", "ind")),
            (("intMinus", "con"), ("intMinus", "non")),
            (("intMinus", "ind"), ("intMinus", "non"))
            ]



add_stat_annotation(g, x="quantile", y="Activity", hue="nodeClass", data=CreUni,  test="Mann-Whitney", box_pairs=boxPairs)

fig.savefig("DirectionalityvsActivity.pdf")


fig = plt.figure(figsize=[6,9])
g = sns.boxplot(x="nodeClass", y="Activity", data=CreUni[CreUni["quantile"] == "biDirectional"], palette=colorPalette)
g.set_yscale("log")
boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]



add_stat_annotation(g, x="nodeClass", y="Activity", data=CreUni[CreUni["quantile"] == "biDirectional"],  test="Mann-Whitney", box_pairs=boxPairs)

fig.savefig("DirectionalityvsActivityBiDirectional.pdf")





sns.scatterplot(x="Activity", y="D.DHT", hue="nodeClass", data=CreUni[CreUni["quantile"] == "biDirectional"])
plt.xscale("log")





# nodeClass = ("tsP", "tsM", "con", "ind", "non")
# for node in nodeClass:
#     sns.distplot(CreUni[CreUni["nodeClass"] == node]["D.DHT"], kde_kws={"label": node}, bins=20)
#     plt.xlim([0,1])







CreUniDir = CrePair[((CrePair["D.DHT"] > 0.9) & (CrePair["M.DHT"] == 1)) | ((CrePair["D.DHT"] < 0.1) & (CrePair["P.DHT"] == 1))]
CreUniwP = CrePair[((CrePair["D.DHT"] > 0.9) & (CrePair["M.DHT"] != 1)) | ((CrePair["D.DHT"] < 0.1) & (CrePair["P.DHT"] != 1))]


CreBiDir = CrePair[(CrePair["D.DHT"] > 0.35) & (CrePair["D.DHT"] < 0.65)]


CreIntLow =  CreUni[((CreUni["D.DHT"] < 0.35) & (CreUni["D.DHT"] > 0.1))]

CreIntHig =  CreUni[((CreUni["D.DHT"] > 0.65) & (CreUni["D.DHT"] < 0.9))]


g = sns.catplot(y="D.DHT", x="nodeClass", data=CreIntLow, kind="violin")
g.fig.get_axes()[0].set_ylim([0,1])
boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non")
        ]
add_stat_annotation(g.fig.get_axes()[0], data=CreIntLow, x="nodeClass", y="D.DHT",  test="Mann-Whitney", box_pairs=boxPairs)


g = sns.catplot(y="D.DHT", x="nodeClass", data=CreIntHig, kind="violin")
g.fig.get_axes()[0].set_ylim([0,1])
boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non")
        ]
add_stat_annotation(g.fig.get_axes()[0], data=CreIntHig, x="nodeClass", y="D.DHT",  test="Mann-Whitney", box_pairs=boxPairs)




nodeClass = ("tss", "con", "ind", "non")
for node in nodeClass:
    sns.distplot(CrePair[CrePair["nodeClass"] == node]["D.DHT"], kde_kws={"label": node}, bins=20)
    plt.xlim([0,1])





# Uni Directionals
g = sns.catplot(y="D.DHT", x="nodeClass", data=CreUniDir, kind="violin")
g = sns.catplot(y="D.DHT", x="nodeClass", data=CreUniwP, kind="violin")


# Uni Directionals
g = sns.catplot(y="D.DHT", x="nodeClass", data=CreBiDir, kind="violin")

g = sns.catplot(y="D.DHT", x="nodeClass", data=CreInt, kind="violin")


g = sns.boxplot(y="D.DHT", x="nodeClass", data=CreInt)





sns.scatterplot(x="P.DHT.Log", y="M.DHT.Log", data=CreBiDir, hue="nodeClass")






CrePlu = CrePair[(CrePair["nodeClass"] == "tsP") | (CrePair["nodeClass"] == "coP" )| (CrePair["nodeClass"] == "inP") | (CrePair["nodeClass"] == "noP")]
CreMin = CrePair[(CrePair["nodeClass"] == "tsM") | (CrePair["nodeClass"] == "coM" )| (CrePair["nodeClass"] == "inM") | (CrePair["nodeClass"] == "noM")]



g = sns.catplot(y="D.DHT", x="nodeClass", data=CreUni, kind="violin")
# g.fig.get_axes()[0].set_yscale('log')

boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non")
        ]


add_stat_annotation(g.fig.get_axes()[0], data=CreUni, x="nodeClass", y="D.DHT",  test="Mann-Whitney", box_pairs=boxPairs)

# Plus

g = sns.catplot(y="D.DHT", x="nodeClass", data=CrePlu, kind="violin")
# g.fig.get_axes()[0].set_yscale('log')

boxPairs = [("tsP", "coP"),
            ("tsP", "inP"),
            ("tsP", "noP")
        ]


add_stat_annotation(g.fig.get_axes()[0], data=CrePlu, x="nodeClass", y="D.DHT",  test="Mann-Whitney", box_pairs=boxPairs)


# Minus

g = sns.catplot(y="D.DHT", x="nodeClass", data=CreMin, kind="violin")
# g.fig.get_axes()[0].set_yscale('log')

boxPairs = [("tsM", "coM"),
            ("tsM", "inM"),
            ("tsM", "noM")
        ]


add_stat_annotation(g.fig.get_axes()[0], data=CreMin, x="nodeClass", y="D.DHT",  test="Mann-Whitney", box_pairs=boxPairs)







#
# u = CrePair["directionValueDHT"].mean()
# SD = CrePair["directionValueDHT"].std()
# CrePair["Z.DHT"] = (CrePair["directionValueDHT"] - u)/SD
#
# u = CrePair["directionValueDMSO"].mean()
# SD = CrePair["directionValueDMSO"].std()
# CrePair["Z.DMSO"] = (CrePair["directionValueDMSO"] - u)/SD
#
# print(CrePair.sort_values(by=['Z.DHT']))

# CrePairDE = CrePair[(abs(CrePair["a.DHT"]) > 1) & (abs(CrePair["b.DHT"]) > 1)]
CrePair["logDHT"] = np.log(CrePair["directionValueDHT"] + 1)
CrePair["absLfc"] = abs(CrePair["Lfc"])

nodeClass = (("tsP", "#f94144"),
             ("tsM", "#577590"),
             ("noM", "#43aa8b"),
             ("noP", "#90be6d"),
             ("coM", "#f8961e"),
             ("coP", "#f3722c"),
             ("inM", "#f8ad9d"),
             ("inP", "#f4978e"))


nodeClass = ("tsP", "tsM", "con", "ind", "non")
for node in nodeClass:
    sns.distplot(CrePair[CrePair["nodeClass"] == node]["Lfc"], kde_kws={"label": node}, bins=20)
    plt.xlim([-5,5])


CrePairtss = CrePair[(CrePair["nodeClass"] == "tsP") | (CrePair["nodeClass"] == "tsM")]
CrePairtss["nodeClass"] = "tss"

CrePair = pd.concat([CrePairtss, CrePair])

CrePairDE = CrePair[(CrePair["absLfc"] > 1)]

g = sns.catplot(y="absLfc", x="nodeClass", data=CrePairDE, kind="violin")
# g.fig.get_axes()[0].set_yscale('log')



boxPairs = [("coP", "tsP"),
            ("inP", "tsP"),
            ("noP", "tsP"),
            ("tsP", "tsM"),
            ("coM", "tsM"),
            ("inM", "tsM"),
            ("noM", "tsM"),
            ("tss", "con"),
            ("tss", "ind"),
            ("tss", "non"),
            ("tss", "tsP"),
            ("tss", "tsM")
        ]


add_stat_annotation(g.fig.get_axes()[0], data=CrePairDE, x="nodeClass", y="absLfc",  test="Mann-Whitney", box_pairs=boxPairs)

plt.savefig("directionViolin.pdf")

























y = "directionValueDHT"
x = "nodeClass"
# y = "BoundTF"
hue = "Strand"
hue_order=['+', '-']

fig = plt.figure(figsize=[9,10])
# colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#FFA96B", "#75E3C1"]
# # colorPalette = ["#FFA96B", "#75E3C1"]
# ax = sns.distplot(CrePair[y], bins=200, fit=stats.gamma)
# ax.set_xlim([-20,20])

nodeClass = ("tsP", "tsM", "non", "con", "ind")
for node in nodeClass:
    sns.distplot(CrePair[CrePair["nodeClass"] == node]["Lfc"], kde_kws={"label": node}, bins=20, hist=False)
    plt.xlim([-6,6])



ax = sns.catplot(y="Lfc", x="nodeClass", data=CrePair, kind="violin")


boxPairs = [("con", "ind"),
            ("con", "non"),
            ("con", "tsP"),
            ("con", "tsM"),
            ("ind", "non"),
            ("ind", "tsP"),
            ("ind", "tsM"),
            ("non", "tsP"),
            ("non", "tsM"),
            ("tsP", "tsM")]

add_stat_annotation(ax, data=CrePair, x="nodeClass", y="Lfc",  test="Mann-Whitney", box_pairs=boxPairs)


g = sns.distplot(CrePair[["directionValueDHT", "nodeClass"]], hue="nodeClass")
g.set(xlim=[-20,20])

plt.savefig("Directional.pdf")

fig = plt.figure(figsize=[9,10])



sns.scatterplot()

ax = sns.scatterplot(x="aLog", y="bLog", data=tips)

y = "b.DHT"
x = "nodeClass"
# y = "BoundTF"
# hue = "Strand"
# hue_order=['+', '-']


colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#FFA96B", "#75E3C1"]
# colorPalette = ["#FFA96B", "#75E3C1"]
ax = sns.boxplot(y=y, x=x, data=CrePair, palette=colorPalette)
ax.set_ylim([0,100])


y = "a.DHT"
x = "nodeClass"
# y = "BoundTF"
# hue = "Strand"
# hue_order=['+', '-']

fig = plt.figure(figsize=[9,10])
colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#FFA96B", "#75E3C1"]
# colorPalette = ["#FFA96B", "#75E3C1"]
ax = sns.boxplot(y=y, x=x, data=CrePair, palette=colorPalette)
ax.set_ylim([0,100])


sns.lmplot(x="aLog", y="bLog", data=CrePair, hue="nodeClass")

sns.lmplot(x="aLog", y="bLog", data=CrePairDE, hue="nodeClass")


fig = plt.figure(figsize=[5,10])
ax=sns.boxplot(y="Lfc", x="nodeClass", data=CrePair)


boxPairs = [("con", "ind"),
            ("con", "non"),
            ("con", "tsP"),
            ("con", "tsM"),
            ("ind", "non"),
            ("ind", "tsP"),
            ("ind", "tsM"),
            ("non", "tsP"),
            ("non", "tsM"),
            ("tsP", "tsM")]

add_stat_annotation(ax, data=CrePair, x="nodeClass", y="Lfc", loc="inside", test="Mann-Whitney", box_pairs=boxPairs)
