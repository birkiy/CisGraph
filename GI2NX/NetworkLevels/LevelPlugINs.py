

from Functions.PlugINs import *


C = pickle.load(open("~/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open("~/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open("~/Data/tmpData/GraphsGData.p", "rb" ))


conBed = {}
readBed(conBed, "~/Data/Regions/cons-arbs.bed")

nonBed = {}
readBed(nonBed, "~/Data/Regions/Non-Active-ARBS.bed")

indBed = {}
readBed(indBed, "~/Data/Regions/ind-arbs.bed")

homBed = {}
readBed(homBed, "~/Data/Regions/promoters_ann_5kb.bed")

tadBed = {}
readBed(tadBed, "~/Data/Regions/TAD.C42B.bed")

chrBed = {}
readBed(chrBed, "~/Data/Regions/hg19.Chr.bed")

beds = {"con": conBed,
        "non": nonBed,
        "ind": indBed,
        "hom": homBed,
        "tad": tadBed,
        "chr": chrBed
        }



colorEdgesPlugIN(G)
colorEdgesPlugIN(T)
colorEdgesPlugIN(C)
print("Colors are added!")


rangePlugIN(G, beds)
rangePlugIN(T, beds)
rangePlugIN(C, beds)
print("Ranges are added!")

fileCSV = "~/Data/DEG/GSE64529_diffexpr-results.csv"

logFCPlugIN(G, fileCSV)
print("LogFC are added!")
