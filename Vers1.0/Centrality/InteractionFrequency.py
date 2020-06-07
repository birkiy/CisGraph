
from Functions.PlotFunctions import *

hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)





e = {}
for edge in G.edges():
    #print(edge)
    edgeType = list(G.edges()[edge]["edgeType"].keys())[0]
    if edgeType not in e.keys():
        e[edgeType] = 1
    else:
        e[edgeType] += 1


edgeN = len(G.edges())
nodeN = len(G.nodes())

D = (2*edgeN)/((nodeN)*(nodeN+1))


edTy = ["upP-upP", "upP-dwP", "upP-con", "upP-ind", "upP-non", "upP-oth",
        "dwP-upP", "dwP-dwP", "dwP-con", "dwP-ind", "dwP-non", "dwP-oth",
        "con-upP", "con-dwP", "con-con", "con-ind", "con-non", "con-oth",
        "ind-upP", "ind-dwP", "ind-con", "ind-ind", "ind-non", "ind-oth",
        "non-upP", "non-dwP", "non-con", "non-ind", "non-non", "non-oth",
        "oth-upP", "oth-dwP", "oth-con", "oth-ind", "oth-non", "oth-oth"]



matIF = []
for edT in edTy:
    c1 = edT.split("-")[0]
    c2 = edT.split("-")[1]
    if c1 == c2:
        V = len([node for node in G.nodes() if G.nodes()[node]["nodeClass"] == c1])
        Emax = (V*(V+1))/2
        if (c1, c2) in e:
            IF = e[(c1, c2)]/ (D*Emax)
        else:
            IF = e[(c2, c1)]/ (D*Emax)
    else:
        V = len([node for node in G.nodes() if G.nodes()[node]["nodeClass"] == c1])
        U = len([node for node in G.nodes() if G.nodes()[node]["nodeClass"] == c2])
        Emax = V*U
        if (c1, c2) in e:
            IF = e[(c1, c2)]/ (D*Emax)
        elif (c2, c1) in e:
            IF = e[(c2, c1)]/ (D*Emax)
        else:
            IF = 0
    matIF += [IF]

matIF = np.array(matIF).reshape(6,6)


fig = plt.figure(figsize=[6,4.5])
ax = sns.heatmap(matIF, cmap=cm,  xticklabels=["upP", "dwP", "con", "ind", "non", "oth"],  yticklabels=["upP", "dwP", "con", "ind", "non", "oth"])
bottom, top = ax.get_ylim()
# ax.set_ylim(bottom + 0.5, top - 0.5)

fig.savefig(f"{figureRoot}/InteractionFrequency.pdf")
