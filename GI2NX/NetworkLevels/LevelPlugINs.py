

colorPalette = {"enh": "#a3c9a8",

                "con": "#FADE89",
                "ind": "#57A4B1",
                "non": "#B0D894",

                "arp": "#F99FA5",
                "ctc": "#9FE2D4",
                "fox": "#C1A2EF",
                "med": "#9ECCF2",
                "enz": "#F7C29C",
                "ari": "#F9E0A1",
                "pol": "#F29ECC",
                "nmy": "#E2B79D",

                "gen": "#ddd8c4",

                "upP": "#95342c",
                "dwP": "#2c6589",

                "tad": "#84b59f",
                "com": "#69a297" ,
                "chr": "#50808e"
}



from Functions.PlugINs import *
from Functions.Helpers import *



C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))



genBed = {}
readBed(genBed, home + "/Data/Regions/promoters_ann_5kb.bed")


enhBed = {}
readBed(enhBed, home + "/Data/Regions/EnhancerDomain.bed")


tadBed = {}
readBed(tadBed, home + "/Data/Regions/TAD.C42B.bed")

chrBed = {}
readBed(chrBed, home + "/Data/Regions/hg19.Chr.bed")





beds = {"gen": genBed,
        "enh": enhBed,
        "tad": tadBed,
        "chr": chrBed,

        }

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



conBed = {}
readBed(conBed, home + "/Data/Regions/cons-arbs.bed")
conBed = {k: v for k,v in sorted(conBed.items(), key=lambda kv: kv[1][1])}

nonBed = {}
readBed(nonBed, home + "/Data/Regions/Non-Active-ARBS.bed")
nonBed = {k: v for k,v in sorted(nonBed.items(), key=lambda kv: kv[1][1])}

indBed = {}
readBed(indBed, home + "/Data/Regions/ind-arbs.bed")
indBed = {k: v for k,v in sorted(indBed.items(), key=lambda kv: kv[1][1])}


Beds = {"con": conBed,
        "ind": indBed,
        "non": nonBed
        }



enhancerFunctionPlugIN(G, Beds, colorPalette)




colorEdgesPlugIN(G, colorPalette)
colorEdgesPlugIN(T, colorPalette)
colorEdgesPlugIN(C, colorPalette)
print("Edges are colored!")


arpBed = {}
readBed(arpBed, home + "/Data/Regions/ProteinBounds/AR.bed")
arpBed = {k: v for k,v in sorted(arpBed.items(), key=lambda kv: kv[1][1])}

ctcBed = {}
readBed(ctcBed, home + "/Data/Regions/ProteinBounds/CTCF.bed")
ctcBed = {k: v for k,v in sorted(ctcBed.items(), key=lambda kv: kv[1][1])}

foxBed = {}
readBed(foxBed, home + "/Data/Regions/ProteinBounds/FOXA1.bed")
foxBed = {k: v for k,v in sorted(foxBed.items(), key=lambda kv: kv[1][1])}

medBed = {}
readBed(medBed, home + "/Data/Regions/ProteinBounds/MED1.bed")
medBed = {k: v for k,v in sorted(medBed.items(), key=lambda kv: kv[1][1])}

enzBed = {}
readBed(enzBed, home + "/Data/Regions/ProteinBounds/EZH2.bed")
enzBed = {k: v for k,v in sorted(enzBed.items(), key=lambda kv: kv[1][1])}

ariBed = {}
readBed(ariBed, home + "/Data/Regions/ProteinBounds/ARID.bed")
ariBed = {k: v for k,v in sorted(ariBed.items(), key=lambda kv: kv[1][1])}

polBed = {}
readBed(polBed, home + "/Data/Regions/ProteinBounds/POLR2.bed")
polBed = {k: v for k,v in sorted(polBed.items(), key=lambda kv: kv[1][1])}


nmyBed = {}
readBed(nmyBed, home + "/Data/Regions/ProteinBounds/NMYC.bed")
nmyBed = {k: v for k,v in sorted(nmyBed.items(), key=lambda kv: kv[1][1])}




Beds = {"arp":arpBed,
        "ctc":ctcBed,
        "fox":foxBed,
        "med":medBed,
        "enz":enzBed,
        "ari":ariBed,
        "pol":polBed,
        "nmy":nmyBed

}

P = nx.Graph()
powerNodePlugIN(G, P, Beds, colorPalette)




pickle.dump(C,open(home + "/Data/tmpData/GraphsCData.p", "wb" ))
pickle.dump(T,open(home + "/Data/tmpData/GraphsTData.p", "wb" ))
pickle.dump(G,open(home + "/Data/tmpData/GraphsGData.p", "wb" ))
pickle.dump(P,open(home + "/Data/tmpData/GraphsPData.p", "wb" ))