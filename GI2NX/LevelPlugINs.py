

from Functions.PlugINs import *


C = pickle.load(open("~/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open("~/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open("~/Data/tmpData/GraphsGData.p", "rb" ))


conBed = {}
readBed(conBed, "~/Data/regions/cons-arbs.bed")

nonBed = {}
readBed(nonBed, "~/Data/regions/Non-Active-ARBS.bed")

indBed = {}
readBed(indBed, "~/Data/regions/ind-arbs.bed")

homBed = {}
readBed(homBed, "~/Data/regions/promoters_ann_5kb.bed")

tadBed = {}
readBed(tadBed, "~/Data/regions/TAD.C42B.bed")

chrBed = {}
readBed(chrBed, "~/Data/regions/hg19.Chr.bed")

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



rangePlugIN(G, beds)
rangePlugIN(T, beds)
rangePlugIN(C, beds)


fileCSV = "~/Data/DEG/GSE64529_diffexpr-results.csv"

logFCPlugIN(G, fileCSV)
