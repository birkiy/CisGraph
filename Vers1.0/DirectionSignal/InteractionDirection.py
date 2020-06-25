
from Functions.Packages import *
from Functions.PlotFunctions import *

G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))

cF = {"+": set(), "-": set()}
cR = {"+": set(), "-": set()}
cD = {"+": set(), "-": set()}
exDF = pd.DataFrame()
for edge in G.edges():
    aEdge = edge[0]
    bEdge = edge[1]
    if G.nodes[aEdge]["lvlP"] == 1 and G.nodes[aEdge]["lvlM"] == 1 and G.nodes[bEdge]["lvlP"] == 1 and G.nodes[bEdge]["lvlM"] == 1:
        continue
    if aEdge == bEdge:
        continue
    aCl = G.nodes[aEdge]["nodeClass"]
    bCl = G.nodes[bEdge]["nodeClass"]
    mAP = G.nodes[aEdge]["lvlP"]
    mAM = G.nodes[aEdge]["lvlM"]
    mBP = G.nodes[bEdge]["lvlP"]
    mBM = G.nodes[bEdge]["lvlM"]
    d = G.edges[(aEdge, bEdge)]["distance"]
    if d != "INF":
        if int(d) < 0:
            dType = "-"
        elif int(d) > 0:
            dType = "+"
    else:
        continue
    if G.nodes[aEdge]["D"] > 0.75:
        dA = "F"
        cF[dType].add(aEdge)
    elif G.nodes[aEdge]["D"] < 0.25:
        dA = "R"
        cR[dType].add(aEdge)
    else:
        dA = "D"
        cD[dType].add(aEdge)
    if G.nodes[bEdge]["D"] > 0.75:
        cF[dType].add(bEdge)
        dB = "F"
    elif G.nodes[bEdge]["D"] < 0.25:
        dB = "R"
        cR[dType].add(bEdge)
    else:
        dB = "D"
        cD[dType].add(bEdge)
    tmpDi = pd.DataFrame({"Distance": dType, "mAP": mAP, "mAM": mAM, "mBP": mBP, "mBM": mBM, "Type": [(dA, dB)], "d": d, "ClassType": [(aCl, bCl)], "nodeA": aEdge, "nodeB": bEdge})
    exDF = pd.concat([exDF, tmpDi])




exDF[['Aclass', 'Bclass']] = pd.DataFrame(exDF['ClassType'].tolist(), index=exDF.index)
exDF = exDF.reset_index(drop=True)


exDF[exDF["Distance"] == "+"].groupby("Type").size()


cl = exDF.groupby("ClassType").size().sort_values()
for c, i in zip(cl, cl.index):
    print(c, i)



Converted = pd.DataFrame()

for i in range(len(exDF)):
    Atmp = exDF.loc[i,"Aclass"]
    Btmp = exDF.loc[i,"Bclass"]
    dtmp = exDF.loc[i,"Distance"]
    if (d == "-"):
        Converted = pd.concat([Converted, exDF.loc[i,:]])
    else:
        typeTmp = exDF.loc[i,"Type"]
        tmpDF = pd.DataFrame({"Distance": "-",
                              "mAP": exDF.loc[i,"mBP"],
                              "mAM": exDF.loc[i,"mBM"],
                              "mBP": exDF.loc[i,"mAP"],
                              "mBM": exDF.loc[i,"mAM"],
                              "Type": [(typeTmp[1], typeTmp[0])],
                              "d": int(exDF.loc[i,"d"])*-1,
                              "ClassType": [(Btmp, Atmp)],
                              "Aclass": Btmp,
                              "Bclass": Atmp})
        Converted = pd.concat([Converted, tmpDF])

Converted = Converted.sort_values("ClassType").reset_index().drop("index", axis=1)

SameClasses = Converted[Converted["Aclass"] == Converted["Bclass"]]
DifferentClasses = Converted[~(Converted["Aclass"] == Converted["Bclass"])]




motor = ["mAP", "mAM", "mBP", "mBM"]

DF0 = pd.DataFrame()
DF0["Tlvl"]= list(Converted[motor[0]])
DF0["Motor"]= motor[0]
DF0["IntType"] = list(Converted["Type"])
DF0["Distance"] = list(Converted["Distance"])
DF0["ClassType"] = list(Converted["ClassType"])

DF1 = pd.DataFrame()
DF1["Tlvl"]= list(Converted[motor[1]])
DF1["Motor"]= motor[1]
DF1["IntType"] = list(Converted["Type"])
DF1["Distance"] = list(Converted["Distance"])
DF1["ClassType"] = list(Converted["ClassType"])

DF2 = pd.DataFrame()
DF2["Tlvl"]= list(Converted[motor[2]])
DF2["Motor"]= motor[2]
DF2["IntType"] = list(Converted["Type"])
DF2["Distance"] = list(Converted["Distance"])
DF2["ClassType"] = list(Converted["ClassType"])

DF3 = pd.DataFrame()
DF3["Tlvl"]= list(Converted[motor[3]])
DF3["Motor"]= motor[3]
DF3["IntType"] = list(Converted["Type"])
DF3["Distance"] = list(Converted["Distance"])
DF3["ClassType"] = list(Converted["ClassType"])


DF = pd.concat([DF0, DF1, DF2, DF3])
DF = DF.sort_values("IntType")


nodeClasses = ["con", "ind", "non", "oth", "pro", "upP", "dwP", "tsP", "uPP", "dPP", "tsM", "uMP", "dMP"]
nodeClasses = ["con", "ind", "non", "oth", "uPP", "dPP", "uMP", "dMP"]
DirTy = [(a, b) for a in nodeClasses for b in nodeClasses]
DF = DF[DF["ClassType"].isin(DirTy)]

DF["logTlvl"] = np.log(DF["Tlvl"] + 1)




hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)





fig = plt.figure(figsize=[20,20])
gs = gridspec.GridSpec(ncols=4, nrows=4)
motors = ["mAP", "mAM", "mBP", "mBM"]
for i, motorA in enumerate(motors):
    for j, motorB in enumerate(motors):
        p = (i,j)
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorA)]["Tlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(8,8)
        cLmeanB = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorB)]["Tlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(8,8)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=nodeClasses,  yticklabels=nodeClasses, center=1)
        plt.title(f"{motorA}/{motorB}")


fig.savefig(f"{figureRoot}/HeatFC.paired.-.pdf")




fig = plt.figure(figsize=[50, 6])

gs = gridspec.GridSpec(ncols=1, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF, x="ClassType", y="logTlvl", hue="Motor", hue_order=["mAP", "mAM", "mBP", "mBM"])
# ax1.set_yscale("log")
plt.title("- aligned")
plt.xticks(rotation=45)

fig.savefig(f"{figureRoot}/ClassTypeMotor.pdf")


























fig = plt.figure(figsize=[15,6])

gs = gridspec.GridSpec(ncols=1, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF, x="Distance", y="logTlvl", hue="Motor", hue_order=["mAP", "mAM", "mBP", "mBM"])
# ax1.set_yscale("log")
plt.title("+")
plt.xticks(rotation=45)

fig.savefig(f"{figureRoot}/DistanceMotor.pdf")

y_vars="logTlvl"
x_vars=[("Distance",(0,0)), ("Motor", (0,1)), ("ClassType", (1,0)), ("IntType", (1,1))]
fig = plt.figure(figsize=[25,25])
gs = gridspec.GridSpec(ncols=2, nrows=2)
for x, p in x_vars:
    fig.add_subplot(gs[p])
    sns.boxplot(data=DF, x=x, y=y_vars)

fig.savefig(f"{figureRoot}/PairPlot.pdf")


hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)

fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP",(0,0)), ("mAM", (0,1)), ("mBP", (1,0)), ("mBM", (1,1))]
for motor, p in motors:
    fig.add_subplot(gs[p])
    cLmean = []
    for ty in DirTy:
        tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motor)]["logTlvl"])
        m = sum(tmpL) / (len(tmpL) +1)
        cLmean += [m]
    cLmean = np.array(cLmean).reshape(8,8)
    ax = sns.heatmap(cLmean, cmap=cm,  xticklabels=nodeClasses,  yticklabels=nodeClasses)
    plt.title(motor)


fig.savefig(f"{figureRoot}/Heat.pdf")




hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)

fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=1)
motors=[("mAP", "mBP"), ("mAM", "mBM")]
for p, tups in enumerate(motors):
    motorA = tups[0]
    motorB = tups[1]
    fig.add_subplot(gs[p])
    cLmeanA = []
    for ty in DirTy:
        tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorA)]["logTlvl"])
        m = (sum(tmpL) +1) / (len(tmpL) +1)
        cLmeanA += [m]
    cLmeanA = np.array(cLmeanA).reshape(8,8)
    cLmeanB = []
    for ty in DirTy:
        tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorB)]["logTlvl"])
        m = (sum(tmpL) +1) / (len(tmpL) +1)
        cLmeanB += [m]
    cLmeanB = np.array(cLmeanB).reshape(8,8)
    cLmean = cLmeanA / cLmeanB
    ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=nodeClasses,  yticklabels=nodeClasses, center=1)
    plt.title(tups)


fig.savefig(f"{figureRoot}/HeatFC.pdf")







fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP", "mBP"), ("mAM", "mBM")]
distances = ["+", "-"]
for i, tups in enumerate(motors):
    for j, d in enumerate(distances):
        p = (i,j)
        motorA = tups[0]
        motorB = tups[1]
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorA) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(8,8)
        cLmeanB = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorB) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(8,8)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=nodeClasses,  yticklabels=nodeClasses, center=1)
        plt.title(f"{tups[0]}/{tups[1]} : {d}")


fig.savefig(f"{figureRoot}/HeatFC.D.pdf")




fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP", "mAM"), ("mBP", "mBM")]
distances = ["+", "-"]
for i, tups in enumerate(motors):
    for j, d in enumerate(distances):
        p = (i,j)
        motorA = tups[0]
        motorB = tups[1]
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorA) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(8,8)
        cLmeanB = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorB) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(8,8)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=nodeClasses,  yticklabels=nodeClasses, center=1)
        plt.title(f"{tups[0]}/{tups[1]} : {d}")


fig.savefig(f"{figureRoot}/HeatFC.D2.pdf")





fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP", "mBM"), ("mAM", "mBP")]
distances = ["+", "-"]
for i, tups in enumerate(motors):
    for j, d in enumerate(distances):
        p = (i,j)
        motorA = tups[0]
        motorB = tups[1]
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorA) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(8,8)
        cLmeanB = []
        for ty in DirTy:
            tmpL = list(DF[(DF["ClassType"] == ty) & (DF["Motor"] == motorB) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(8,8)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=nodeClasses,  yticklabels=nodeClasses, center=1)
        plt.title(f"{tups[0]}/{tups[1]} : {d}")


fig.savefig(f"{figureRoot}/HeatFC.D3.pdf")






fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP", "mBM"), ("mAM", "mBP")]

IntTypes = [(a, b) for a in ["F", "R", "D"] for b in ["F", "R", "D"]]
distances = ["+", "-"]
for i, tups in enumerate(motors):
    for j, d in enumerate(distances):
        p = (i,j)
        motorA = tups[0]
        motorB = tups[1]
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in IntTypes:
            tmpL = list(DF[(DF["IntType"] == ty) & (DF["Motor"] == motorA) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(3,3)
        cLmeanB = []
        for ty in IntTypes:
            tmpL = list(DF[(DF["IntType"] == ty) & (DF["Motor"] == motorB) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(3,3)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=["F", "R", "D"],  yticklabels=["F", "R", "D"], center=1)
        plt.title(f"{tups[0]}/{tups[1]} : {d}")


fig.savefig(f"{figureRoot}/HeatFC.D4.pdf")





fig = plt.figure(figsize=[20,12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
motors=[("mAP", "mBP"), ("mAM", "mBM")]

IntTypes = [(a, b) for a in ["F", "R", "D"] for b in ["F", "R", "D"]]
distances = ["+", "-"]
for i, tups in enumerate(motors):
    for j, d in enumerate(distances):
        p = (i,j)
        motorA = tups[0]
        motorB = tups[1]
        fig.add_subplot(gs[p])
        cLmeanA = []
        for ty in IntTypes:
            tmpL = list(DF[(DF["IntType"] == ty) & (DF["Motor"] == motorA) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanA += [m]
        cLmeanA = np.array(cLmeanA).reshape(3,3)
        cLmeanB = []
        for ty in IntTypes:
            tmpL = list(DF[(DF["IntType"] == ty) & (DF["Motor"] == motorB) & (DF["Distance"] == d)]["logTlvl"])
            m = (sum(tmpL) +1) / (len(tmpL) +1)
            cLmeanB += [m]
        cLmeanB = np.array(cLmeanB).reshape(3,3)
        cLmean = cLmeanA / cLmeanB
        ax = sns.heatmap(cLmean, cmap=sns.color_palette("coolwarm", 200),  xticklabels=["F", "R", "D"],  yticklabels=["F", "R", "D"], center=1)
        plt.title(f"{tups[0]}/{tups[1]} : {d}")


fig.savefig(f"{figureRoot}/HeatFC.D5.pdf")







fig = plt.figure(figsize=[200,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF[DF["Distance"] == "+"], x="ClassType", y="Tlvl", hue="Motor", order=DirTy, hue_order=["mAP", "mAM", "mBP", "mBM"])
ax1.set_yscale("log")
plt.title("+")
plt.xticks(rotation=45)










fig = plt.figure(figsize=[200,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF[DF["Distance"] == "+"], x="ClassType", y="Tlvl", hue="Motor", order=DirTy, hue_order=["mAP", "mAM", "mBP", "mBM"])
ax1.set_yscale("log")
plt.title("+")
plt.xticks(rotation=45)

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=DF[DF["Distance"] == "-"], x="ClassType", y="Tlvl", hue="Motor", order=DirTy, hue_order=["mAP", "mAM", "mBP", "mBM"])
ax2.set_yscale("log")
plt.title("-")
plt.xticks(rotation=45)

fig.savefig(f"{figureRoot}/ClassTypeMotor.pdf")



DirTy2 = [(a, b      Distance       mAP       mAM       mBP        mBM    Type d    ClassType Aclass Bclass
23           -  126.6860   42.3402    1.0000     1.0000  (F, F)     (con, ind)    con    ind
24           -   22.0872    1.0000    1.0000     1.0000  (D, D)     (con, ind)    con    ind
25           -  142.3270   14.7801    1.0000     1.0000  (F, D)     (con, ind)    con    ind
26           -   63.9824   42.3403   42.8952    56.3959  (F, D)     (con, ind)    con    ind
27           -    1.0000  198.7440    1.0000     1.0000  (R, D)     (con, ind)    con    ind
...        ...       ...       ...       ...        ...     ... ..         ...    ...    ...
22153        -    0.0000    0.0000    1.0000     1.0000  (D, D)     (upP, uMP)) for a in ["F", "R", "D"] for b in ["F", "R", "D"]]

DF = DF.sort_values("ClassType")

fig = plt.figure(figsize=[300,20])

gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF[DF["Distance"] == "+"], x="IntType", y="Tlvl", hue="ClassType", order=DirTy2, hue_order=DirTy)
ax1.set_yscale("log")
plt.title("+")
plt.xticks(rotation=45)

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=DF[DF["Distance"] == "-"], x="IntType", y="Tlvl", hue="ClassType", order=DirTy2, hue_order=DirTy)
ax2.set_yscale("log")
plt.title("-")
plt.xticks(rotation=45)

fig.savefig(f"{figureRoot}/IntTypeClassType.pdf")





fig = plt.figure(figsize=[15,15])

gs = gridspec.GridSpec(ncols=2, nrows=2)
fig.add_subplot(gs[0,0])
ax1 = sns.boxplot(data=DF[DF["Motor"] == motor[0]], x="IntType", y="Tlvl", hue="Distance", order=DirTy)
ax1.set_yscale("log")
plt.title(motor[0])

fig.add_subplot(gs[0,1])
ax2 = sns.boxplot(data=DF[DF["Motor"] == motor[0]], x="IntType", y="Tlvl", hue="Distance", order=DirTy)
ax2.set_yscale("log")
plt.title(motor[1])


fig.add_subplot(gs[1,0])
ax1 = sns.boxplot(data=DF[DF["Motor"] == motor[2]], x="IntType", y="Tlvl", hue="Distance", order=DirTy)
ax1.set_yscale("log")
plt.title(motor[2])

fig.add_subplot(gs[1,1])
ax2 = sns.boxplot(data=DF[DF["Motor"] == motor[3]], x="IntType", y="Tlvl", hue="Distance", order=DirTy)
ax2.set_yscale("log")
plt.title(motor[3])


fig.savefig(f"{figureRoot}/mapDist.pdf")
