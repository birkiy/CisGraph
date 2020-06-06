


from Functions.DeltaCon import *



ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

d, w, E, G1, G2= NX2deltaCon(ethG, dhtG, returnG=True, e=10**-5)



Impact = pd.DataFrame()
Impact["w"] = dict(G1.nodes(data="w")).values()
Impact["nodeClass"] = dict(G1.nodes(data="nodeClass")).values()
Impact["node"] = dict(G1.nodes(data="nodeClass")).keys()


Arbs = Impact[~((Impact["nodeClass"] == "upP") | (Impact["nodeClass"] == "dwP"))]

boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "oth"),
            ("ind", "non"),
            ("ind", "oth"),
            ("non", "oth"))


fig = plt.figure(figsize=(5,8))
ax = sns.boxplot(x="nodeClass", y="w", data=Arbs, palette=colorPalette, width = 0.5)
add_stat_annotation(ax,x="nodeClass", y="w", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
plt.title("DHT vs. EtOH Co-accesibility\n (DeltaCon Scores)")
fig.savefig(f"{figureRoot}/deltaConE5.DHTvsEtOH.pdf")


plt.close("all")
