
from Functions.Utils import *


dhtG = pickle.load(open(f"{dataRoot}/PickleData/GraphsG.DHT.Data.p", "rb" ))
ethG = pickle.load(open(f"{dataRoot}/PickleData/GraphsG.EtOH.Data.p", "rb" ))


Dir = {"PA": {"+": 0, "-": 0},"MA": {"+": 0, "-": 0}, "OA": {"+": 0, "-": 0}, "PB": {"+": 0, "-": 0},"MB": {"+": 0, "-": 0}, "OB": {"+": 0, "-": 0}}
Dir = {("P", "M"): {"+": 0, "-": 0}, ("M", "P"): {"+": 0, "-": 0},
       ("M", "O"): {"+": 0, "-": 0}, ("O", "P"): {"+": 0, "-": 0},
       ("O", "M"): {"+": 0, "-": 0}, ("P", "O"): {"+": 0, "-": 0},
       ("O", "O"): {"+": 0, "-": 0}, ("P", "P"): {"+": 0, "-": 0}, ("M", "M"): {"+": 0, "-": 0}}
visited = []
dC = {"+": {"P": 0, "M": 0, "O": 0}, "-": {"P": 0, "M": 0, "O": 0}}
exD = {"+": {("P", "M"): [], ("M", "P"): [],
             ("M", "O"): [], ("O", "P"): [],
             ("O", "M"): [], ("P", "O"): [],
             ("O", "O"): [], ("P", "P"): [], ("M", "M"): []},
       "-": {("P", "M"): [], ("M", "P"): [],
             ("M", "O"): [], ("O", "P"): [],
             ("O", "M"): [], ("P", "O"): [],
             ("O", "O"): [], ("P", "P"): [], ("M", "M"): []}}
exDF = pd.DataFrame()
for edge in dhtG.edges():
    aEdge = edge[0]
    bEdge = edge[1]
    if dhtG.nodes[aEdge]["lvlP"] == 1 and dhtG.nodes[aEdge]["lvlM"] == 1 and dhtG.nodes[bEdge]["lvlP"] == 1 and dhtG.nodes[bEdge]["lvlM"] == 1:
        continue
    if dhtG.nodes[aEdge]["D"] > 0.75:
        dCA = "P"
        exA = dhtG.nodes[aEdge]["lvlP"]
    elif dhtG.nodes[aEdge]["D"] < 0.25:
        dCA = "M"
        exA = dhtG.nodes[aEdge]["lvlM"]
    else:
        dCA = "O"
        exA = dhtG.nodes[aEdge]["lvlP"] - dhtG.nodes[aEdge]["lvlM"]
    if dhtG.nodes[bEdge]["D"] > 0.75:
        dCB = "P"
        exB = dhtG.nodes[bEdge]["lvlP"]
    elif dhtG.nodes[bEdge]["D"] < 0.25:
        dCB = "M"
        exB = dhtG.nodes[bEdge]["lvlM"]
    else:
        dCB = "O"
        exB = dhtG.nodes[aEdge]["lvlP"] - dhtG.nodes[bEdge]["lvlM"]
    d = dhtG[bEdge][aEdge]["distance"]
    if d != "INF":
        if int(d) < 0:
            tmpDi = pd.DataFrame({"Distance": "-", "exA": exA, "exB": exB, "Type": [(dCA, dCB)], "d": d})
            Dir[(dCA, dCB)]["-"] += 1
            dC["-"][dCA] += 1
            dC["-"][dCB] += 1
            exD["-"][(dCA, dCB)] += [exA]
        elif int(d) > 0:
            tmpDi = pd.DataFrame({"Distance": "+", "exA": exA, "exB": exB, "Type": [(dCA, dCB)], "d": d})
            Dir[(dCA, dCB)]["+"] += 1
            dC["+"][dCA] += 1
            dC["+"][dCB] += 1
            exD["+"][(dCA, dCB)] += [exA]
        exDF = pd.concat([exDF, tmpDi])

sumP = 0
sumM = 0
for ty in Dir.keys():
    sumP += Dir[ty]["+"]
    sumM += Dir[ty]["-"]




cF = {"+": set(), "-": set()}
cR = {"+": set(), "-": set()}
cD = {"+": set(), "-": set()}
exDF = pd.DataFrame()
for edge in dhtG.edges():
    aEdge = edge[0]
    bEdge = edge[1]
    if dhtG.nodes[aEdge]["lvlP"] == 1 and dhtG.nodes[aEdge]["lvlM"] == 1 and dhtG.nodes[bEdge]["lvlP"] == 1 and dhtG.nodes[bEdge]["lvlM"] == 1:
        continue
    if aEdge == bEdge:
        continue
    mAP = dhtG.nodes[aEdge]["lvlP"]
    mAM = dhtG.nodes[aEdge]["lvlM"]
    mBP = dhtG.nodes[bEdge]["lvlP"]
    mBM = dhtG.nodes[bEdge]["lvlM"]
    d = dhtG.edges[(aEdge, bEdge)]["distance"]
    if d != "INF":
        if int(d) < 0:
            dType = "-"
        elif int(d) > 0:
            dType = "+"
    else:
        print("hey")
        continue
    if dhtG.nodes[aEdge]["D"] > 0.75:
        dA = "F"
        cF[dType].add(aEdge)
    elif dhtG.nodes[aEdge]["D"] < 0.25:
        dA = "R"
        cR[dType].add(aEdge)
    else:
        dA = "D"
        cD[dType].add(aEdge)
    if dhtG.nodes[bEdge]["D"] > 0.75:
        cF[dType].add(bEdge)
        dB = "F"
    elif dhtG.nodes[bEdge]["D"] < 0.25:
        dB = "R"
        cR[dType].add(bEdge)
    else:
        dB = "D"
        cD[dType].add(bEdge)
    tmpDi = pd.DataFrame({"Distance": dType, "mAP": mAP, "mAM": mAM, "mBP": mBP, "mBM": mBM, "Type": [(dA, dB)], "d": d})
    exDF = pd.concat([exDF, tmpDi])


exDF[exDF["Distance"] == "+"].groupby("Type").size()


motor = ["mAP", "mAM", "mBP", "mBM"]

DF0 = pd.DataFrame()
DF0["Tlvl"]= list(exDF[motor[0]])
DF0["Motor"]= motor[0]
DF0["IntType"] = list(exDF["Type"])
DF0["Distance"] = list(exDF["Distance"])
DF1 = pd.DataFrame()
DF1["Tlvl"]= list(exDF[motor[1]])
DF1["Motor"]= motor[1]
DF1["IntType"] = list(exDF["Type"])
DF1["Distance"] = list(exDF["Distance"])
DF2 = pd.DataFrame()
DF2["Tlvl"]= list(exDF[motor[2]])
DF2["Motor"]= motor[2]
DF2["IntType"] = list(exDF["Type"])
DF2["Distance"] = list(exDF["Distance"])
DF3 = pd.DataFrame()
DF3["Tlvl"]= list(exDF[motor[3]])
DF3["Motor"]= motor[3]
DF3["IntType"] = list(exDF["Type"])
DF3["Distance"] = list(exDF["Distance"])

DF = pd.concat([DF0, DF1, DF2, DF3])
DF = DF.sort_values("Type")

fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF[DF["Distance"] == "+"], x="IntType", y="Tlvl", hue="Motor", order=DirTy)
ax1.set_yscale("log")
plt.title("+")

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=DF[DF["Distance"] == "-"], x="IntType", y="Tlvl", hue="Motor", order=DirTy)
ax2.set_yscale("log")
plt.title("-")

fig.savefig(f"{figureRoot}/map.pdf")





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











DF.groupby("IntType").size() // 4






Emax = (V*(V+1))

Emax = V*U*2


from Functions.PlotFunctions import *


c = {"F": cF, "R": cR, "D": cD}
DirTy = [(a, b) for a in ["F", "R", "D"] for b in ["F", "R", "D"]]
Dir = exDF.groupby(["Type", "Distance"]).size().reset_index()
matP = []
matM = []
for ty in DirTy:
    if ty[0] == ty[1]:
        VP = len(c[ty[0]]["+"])
        EmaxP = (VP*(VP+1))/2
        VM = len(c[ty[0]]["-"])
        EmaxM = (VM*(VM+1))/2
    else:
        VP = len(c[ty[0]]["+"])
        UP = len(c[ty[1]]["+"])
        EmaxP = VP*UP
        VM = len(c[ty[0]]["-"])
        UM = len(c[ty[1]]["-"])
        EmaxM = VM*UM
    p = int(Dir[(Dir["Type"] == ty) & (Dir["Distance"] == "+")].iloc[:,2]) / EmaxP
    m = int(Dir[(Dir["Type"] == ty) & (Dir["Distance"] == "-")].iloc[:,2]) / EmaxM
    matP += [p]
    matM += [m]

matP = np.array(matP).reshape(3,3)
matM = np.array(matM).reshape(3,3)



hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)

fig = plt.figure(figsize=[9,5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax = sns.heatmap(matP, cmap=cm,  xticklabels=["F", "R", "D"],  yticklabels=["F", "R", "D"])
# bottom, top = ax.get_ylim()
# ax.set_ylim(bottom + 0.5, top - 0.5)


#
#
# hc = ["#ffbaba", "#ff7b7b", "#ff5252", "#ff0000", "#A70000"]
#
# th = [0, 0.15, 0.5, 0.85, 1]
# cdict = NonLinCdict(th, hc)
# cm = colors.LinearSegmentedColormap('test', cdict)

fig.add_subplot(gs[1])
ax = sns.heatmap(matM, cmap=cm,  xticklabels=["F", "R", "D"],  yticklabels=["F", "R", "D"])

fig.savefig(f"{figureRoot}/DirectionHeat.pdf")
