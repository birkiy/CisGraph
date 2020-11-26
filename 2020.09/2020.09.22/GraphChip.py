


from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))

chipAll = pickle.dump(open(f"{dataRoot}/ChipAll.p", "rb" ))

chipAll.loc[chipAll.index[chipAll["Unnamed: 0"] == "peakno_1466-ARBS"], :]

splitChip = chipAll.set_index("Unnamed: 0").drop(columns=["nodeClass"]).to_dict('split')

corDF = pd.DataFrame()
for i, edge in enumerate(G.edges()):
    aEdge = edge[0]
    bEdge = edge[1]
    if aEdge == bEdge:
        continue
    if aEdge not in splitChip["index"] or bEdge not in splitChip["index"]:
        continue
    aCl = G.nodes[aEdge]["nodeClass"]
    bCl = G.nodes[bEdge]["nodeClass"]
    if aCl in ("tsP", "tsM", "uMP", "uPP", "dPP", "dMP"):
        aCl = "tss"
    if bCl in ("tsP", "tsM", "uMP", "uPP", "dPP", "dMP"):
        bCl = "tss"
    if aCl == "oth" or bCl == "oth":
        continue
    aIdx = splitChip["index"].index(aEdge)
    bIdx = splitChip["index"].index(bEdge)
    aTFs = splitChip["data"][aIdx]
    bTFs = splitChip["data"][bIdx]
    cor, p = scipy.stats.spearmanr(aTFs, bTFs)
    tmpDi = pd.DataFrame({"Edge": [(aEdge, bEdge)], "Cor": cor, "pVal": p, "ClassType": [(aCl, bCl)]})
    corDF = pd.concat([corDF, tmpDi])
    if i % 1000 == 0:
        print(i)

def fdr(p_vals):
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

corDF["qVal"] = fdr(corDF["pVal"])

corDF = corDF[corDF["qVal"] < 0.05]

corDF.["nodeClass"]


fig = plt.figure(figsize=[5,5])
gs = gridspec.GridSpec(ncols=1, nrows=1)
plt.subplots_adjust(bottom=0.4)

ax = sns.boxplot(data=corDF, x="ClassType", y="Cor")
plt.xticks(rotation=45)
fig.savefig(f"{figureRoot}/ChipGraph.pdf")
