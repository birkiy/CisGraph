

from Functions.PlugINs import *
from Functions.Helpers import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))


conBed = {}
readBed(conBed, home + "/Data/Regions/cons-arbs.bed")

nonBed = {}
readBed(nonBed, home + "/Data/Regions/Non-Active-ARBS.bed")

indBed = {}
readBed(indBed, home + "/Data/Regions/ind-arbs.bed")

homBed = {}
readBed(homBed, home + "/Data/Regions/promoters_ann_5kb.bed")

tadBed = {}
readBed(tadBed, home + "/Data/Regions/TAD.C42B.bed")

chrBed = {}
readBed(chrBed, home + "/Data/Regions/hg19.Chr.bed")

beds = {"con": conBed,
        "non": nonBed,
        "ind": indBed,
        "hom": homBed,
        "tad": tadBed,
        "chr": chrBed
        }

colorPalette ={"hom": "#888888",
               "con": "#FADE89",
               "ind": "#57A4B1",
               "non": "#B0D894",
               "tad": "#4c5172",
               "chr": "#5a4c72"
               }


colorEdgesPlugIN(G, colorPalette)
colorEdgesPlugIN(T, colorPalette)
colorEdgesPlugIN(C, colorPalette)
print("Colors are added!")


rangePlugIN(G, beds)
rangePlugIN(T, beds)
rangePlugIN(C, beds)
print("Ranges are added!")


fileCSV = home + "/Data/DEG/GSE64529_diffexpr-results.csv"

logFCPlugIN(G, fileCSV)
print("LogFC are added!")



upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "dwP"]

print(len(upP), len(dwP))



pickle.dump(C,open(home + "/Data/tmpData/GraphsCData.p", "wb" ))
pickle.dump(T,open(home + "/Data/tmpData/GraphsTData.p", "wb" ))
pickle.dump(G,open(home + "/Data/tmpData/GraphsGData.p", "wb" ))
