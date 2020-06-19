
from Functions.PlugINs import *
from Functions.Helpers import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))

ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

proBed = {}
readBed(proBed, f"{dataRoot}/Regions/promoters_ann_5kb.bed")

conBed = {}
readBed(conBed, f"{dataRoot}/Regions/cons-arbs.bed")

indBed = {}
readBed(indBed, f"{dataRoot}/Regions/ind-arbs.bed")

nonBed = {}
readBed(nonBed, f"{dataRoot}/Regions/Non-Active-ARBS.bed")

othBed = {}
readBed(othBed, f"{dataRoot}/Regions/otherARBS.bed")


tsPBed = {}
readBed(tsPBed, "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.4.+.bed")
tsPBed = dict(pd.DataFrame(tsPBed))

tsMBed = {}
readBed(tsMBed, "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.4.-.bed")

G.remove_node("")

beds = {"pro": proBed,
        "con": conBed,
        "ind": indBed,
        "non": nonBed,
        "oth": othBed,
        "tsP": tsPBed,
        "tsM": tsMBed
        }
rangePlugIN(G, beds)

rangePlugIN(ethG, beds)
rangePlugIN(dhtG, beds)
print("Ranges are added!")



GroDF = pickle.load(open(f"{dataRoot}/tmpData/Gro.DF.p", "rb" ))
for node in G.nodes():
    if G.nodes[node]["nodeClass"] == "pro":
        continue
    G.nodes[node]["lvlM"] = GroDF.loc[node, "dmso.-"]
    G.nodes[node]["lvlP"] = GroDF.loc[node, "dmso.+"]
    #
    G.nodes[node]["D"] = GroDF.loc[node, "Directionality"]

print("Directionals are added!")


print("\n")

fileCSV = f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv"
logFCPlugIN(G, fileCSV, colorPalette={"upP": "#000000", "dwP": "#000000"} )
logFCPlugIN(G, fileCSV, colorPalette={"uPP": "#000000", "dPP": "#000000"}, geneClass="tsP")
logFCPlugIN(G, fileCSV, colorPalette={"uMP": "#000000", "dMP": "#000000"}, geneClass="tsM")

logFCPlugIN(ethG, fileCSV, colorPalette={"upP": "#000000", "dwP": "#000000"} )
logFCPlugIN(dhtG, fileCSV, colorPalette={"upP": "#000000", "dwP": "#000000"} )
print("LogFC are added!")


nodesG = list(G.nodes())
for node in nodesG:
    if G.nodes()[node]["nodeClass"] == "pro":
        G.remove_node(node)


nodesG = list(ethG.nodes())
for node in nodesG:
    if ethG.nodes()[node]["nodeClass"] == "pro":
        ethG.remove_node(node)


nodesG = list(dhtG.nodes())
for node in nodesG:
    if dhtG.nodes()[node]["nodeClass"] == "pro":
        dhtG.remove_node(node)




pickle.dump(G,open(f"{dataRoot}/tmpData/GraphsG.Data.p", "wb" ))



print("\n")
print("You have new GENE classes \"upP\" and \"dwP\". \nTheir numbers relatively:")
upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "dwP"]
print(len(upP), len(dwP))



print("\n")
print("You have new GENE classes \"upP\" and \"dwP\". \nTheir numbers relatively:")
upP = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "dwP"]
print(len(upP), len(dwP))

print("\n")
print("You have new GENE classes \"upP\" and \"dwP\". \nTheir numbers relatively:")
upP = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "dwP"]
print(len(upP), len(dwP))


print("You have new ENHANCER classes \"con\" , \"ind\" and \"non\". \nTheir numbers relatively:")
con = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "con"]
ind = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "ind"]
non = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "non"]
oth = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "oth"]
print(len(con), len(ind), len(non), len(oth))
print("\n")


print("You have new ENHANCER classes \"con\" , \"ind\" and \"non\". \nTheir numbers relatively:")
con = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "con"]
ind = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "ind"]
non = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "non"]
oth = [_[0] for _ in ethG.nodes(data="nodeClass") if _[1] == "oth"]
print(len(con), len(ind), len(non), len(oth))
print("\n")


print("You have new ENHANCER classes \"con\" , \"ind\" and \"non\". \nTheir numbers relatively:")
con = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "con"]
ind = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "ind"]
non = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "non"]
oth = [_[0] for _ in dhtG.nodes(data="nodeClass") if _[1] == "oth"]
print(len(con), len(ind), len(non), len(oth))
print("\n")

pickle.dump(G,open(f"{dataRoot}/tmpData/GraphsGData.p", "wb" ))
pickle.dump(ethG,open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "wb" ))
pickle.dump(dhtG,open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "wb" ))
