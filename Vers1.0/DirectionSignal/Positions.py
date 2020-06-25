


from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


distanceDF = pd.DataFrame()
for i, edge in enumerate(G.edges()):
    aEdge = edge[0]
    bEdge = edge[1]
    if aEdge == bEdge:
        continue
    aCl = G.nodes[aEdge]["nodeClass"]
    bCl = G.nodes[bEdge]["nodeClass"]
    d = G.edges[(aEdge, bEdge)]["distance"]
    if d != "INF":
        if int(d) < 0:
            dType = "-"
        elif int(d) > 0:
            dType = "+"
    else:
        continue
    tmpDi = pd.DataFrame({"Position": dType, "Distance": int(d), "ClassType": [(aCl, bCl)]})
    distanceDF = pd.concat([distanceDF, tmpDi])
    if i % 1000 == 0:
        print(i)



distanceDF[['Aclass', 'Bclass']] = pd.DataFrame(distanceDF['ClassType'].tolist(), index=distanceDF.index)
distanceDF = distanceDF.sort_values("Distance").reset_index().drop("index", axis=1)



nodeClasses = ["con", "ind","non","upP","dwP"]

nodeClasses = ["con", "ind", "non", "oth", "uPP", "uMP", "upP"]



Arbs = distanceDF[distanceDF["Aclass"].isin(nodeClasses) & distanceDF["Bclass"].isin(nodeClasses)].reset_index(drop=True)




Arbs.loc[Arbs["Position"] == "-", "Distance"] = Arbs.loc[Arbs["Position"] == "-", "Distance"] * -1
fig = plt.figure(figsize=[10,10])

ax = sns.boxplot(data=Arbs, x="ClassType", y="Distance", hue="Position")
ax.set_yscale("log")
plt.xticks(rotation=45)
fig.savefig(f"{figureRoot}/DistancePosition.pdf")


nodeClasses = [((0,0),"con"),((0,1), "ind"),((0,2),"non"), ((1,0),"upP"),((1,1),"dwP")]

p2m = lambda b: "-" if b == "+" else "+"
t2C = lambda tup: (tup[1], tup[0])
fig = plt.figure(figsize=[15,10])
gs = gridspec.GridSpec(ncols=3, nrows=2)
for gsp, nodeA in nodeClasses:
    Aal = Arbs[Arbs["Aclass"] == nodeA].reset_index(drop=True)
    Bal = Arbs[(Arbs["Bclass"] == nodeA) & (Arbs["Aclass"] != nodeA)].reset_index(drop=True)
    for i in range(len(Bal)):
        p = p2m(Bal.loc[i,"Position"])
        d = Bal.loc[i,"Distance"] * -1
        c = t2C(Bal.loc[i,"ClassType"])
        a = Bal.loc[i,"Bclass"]
        b = Bal.loc[i,"Aclass"]
        tmp = pd.DataFrame({"Position": p, "Distance": int(d), "ClassType": [c], "Aclass": a, "Bclass": b})
        Aal = pd.concat([Aal, tmp])
    Aal = Aal.sort_values("Distance").reset_index(drop=True)
    Aal["logDistance"] = Aal.apply(lambda b: np.log(b["Distance"]+1) if b["Distance"] > 0 else (np.log((b["Distance"]-1)*-1))*-1, axis=1)
    kwargs = {'cumulative': True}
    fig.add_subplot(gs[gsp])
    for i, node in nodeClasses:
        x = Aal[(Aal["Bclass"] == node) ]["logDistance"]
        ax1 = sns.distplot(x, hist=False, kde_kws=kwargs, label=node + ": " + str(len(x)) )
    plt.legend()
    plt.title(nodeA)



fig.savefig(f"{figureRoot}/CDF.Position.pdf")

plt.close("all")


ArbsDif = Arbs.query("Aclass != Bclass" )

p2m = lambda b: "-" if b == "+" else "+"
t2C = lambda tup: (tup[1], tup[0])
fig = plt.figure(figsize=[10,30])
gs = gridspec.GridSpec(ncols=2, nrows=len(nodeClasses))
plt.subplots_adjust(top=0.95, right=0.95, wspace=0.2, hspace=0.3)
for idx, nodeA in enumerate(nodeClasses):
    Aal = ArbsDif[ArbsDif["Aclass"] == nodeA].reset_index(drop=True)
    Bal = ArbsDif[(ArbsDif["Bclass"] == nodeA) & (ArbsDif["Aclass"] != nodeA)].reset_index(drop=True)
    for i in range(len(Bal)):
        p = p2m(Bal.loc[i,"Position"])
        d = Bal.loc[i,"Distance"] * -1
        c = t2C(Bal.loc[i,"ClassType"])
        a = Bal.loc[i,"Bclass"]
        b = Bal.loc[i,"Aclass"]
        tmp = pd.DataFrame({"Position": p, "Distance": int(d), "ClassType": [c], "Aclass": a, "Bclass": b})
        Aal = pd.concat([Aal, tmp])
    Aal = Aal.sort_values("Distance").reset_index(drop=True)
    Aal["logDistance"] = Aal.apply(lambda b: np.log(b["Distance"]+1) if b["Distance"] > 0 else (np.log((b["Distance"]-1)*-1))*-1, axis=1)
    kwargs = {'cumulative': True}
    fig.add_subplot(gs[(idx, 0)])
    for node in nodeClasses:
        x = - Aal[(Aal["Bclass"] == node) ]["logDistance"]
        ax1 = sns.distplot(x, hist=False, kde_kws=kwargs, label=node + ": " + str(len(x)) )
    plt.legend()
    plt.title(nodeA + "CDF")
    fig.add_subplot(gs[(idx), 1])
    for node in nodeClasses:
        x = - Aal[(Aal["Bclass"] == node) ]["Distance"]
        ax2 = sns.distplot(x, hist=False, label=node + ": " + str(len(x)) )
    plt.legend()
    plt.title(nodeA + "PDF")


fig.savefig(f"{figureRoot}/CDFvsPDF.Dif.all.Position.pdf")

plt.close("all")

fig = plt.figure(figsize=[9,4.5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(top=0.95, right=0.95, wspace=0.2, hspace=0.3)

ArbsSam = Arbs.query("Aclass == Bclass").reset_index(drop=True)
ArbsSam["logDistance"] = ArbsSam.apply(lambda b: np.log(b["Distance"]+1) if b["Distance"] > 0 else (np.log((b["Distance"]-1)*-1))*-1, axis=1)

kwargs = {'cumulative': True}
fig.add_subplot(gs[0])
for node in nodeClasses:
    x = abs(ArbsSam[(ArbsSam["Bclass"] == node) ]["logDistance"])
    ax1 = sns.distplot(x, hist=False, kde_kws=kwargs, label=node + ": " + str(len(x)) )

plt.legend()
plt.title("Same Class Interactions." + "CDF")
fig.add_subplot(gs[1])
for node in nodeClasses:
    x = abs(ArbsSam[(ArbsSam["Bclass"] == node) ]["logDistance"])
    ax2 = sns.distplot(x, hist=False, label=node + ": " + str(len(x)) )

plt.legend()
plt.title("Same Class Interactions." + "PDF")


fig.savefig(f"{figureRoot}/CDFvsPDF.Sam.Position.pdf")








D = pd.DataFrame()
for i in range(len(ArbsSam)):
    d = ArbsSam.loc[i, "Distance"]
    c = ArbsSam.loc[i, "Aclass"]
    tmp = pd.DataFrame({"Distance": d, "class": c}, index=[0])
    D = pd.concat([D, tmp])

D = D.sort_values("Distance").reset_index(drop=True)
