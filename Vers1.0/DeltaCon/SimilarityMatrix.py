


from Functions.DeltaCon import *




G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))

ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

d, w, E, G1, G2= NX2deltaCon(ethG, dhtG, returnG=True, e=10**-8)

d1, w1, E1, G11, G21= NX2deltaCon(G, dhtG, returnG=True, e=10**-8)

Impact = pd.DataFrame()
Impact["w"] = dict(G11.nodes(data="w")).values()
Impact["nodeClass"] = dict(G11.nodes(data="nodeClass")).values()
Impact["node"] = dict(G11.nodes(data="nodeClass")).keys()



colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "upP": "#000000",
                "dwP": "#000000"}

boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "dwP"),
            ("con", "upP"),
            ("con", "oth"),
            ("ind", "non"),
            ("ind", "upP"),
            ("ind", "dwP"),
            ("ind", "oth"),
            ("non", "upP"),
            ("non", "oth"),
            ("non", "dwP"),
            ("upP", "dwP"),
            ("upP", "oth"),
            ("oth", "dwP"))



boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "pro"),
            ("con", "oth"),
            ("ind", "non"),
            ("ind", "pro"),
            ("ind", "oth"),
            ("non", "pro"),
            ("non", "oth"),
            ("pro", "oth"),
            ("oth", "pro"))

fig = plt.figure(figsize=(9,12))
ax = sns.violinplot(x="nodeClass", y="w", data=Impact, palette=colorPalette)
add_stat_annotation(ax,x="nodeClass", y="w", data=Impact,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
plt.title("DHT Co-accesibility vs. G\n (DeltacCon Scores)")
fig.savefig(f"{figureRoot}/deltaConE8.GvsDHTwOth.pdf")


plt.close("all")
