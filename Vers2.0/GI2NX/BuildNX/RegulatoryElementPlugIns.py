
from Functions.FeaturePlugIns import *


dhtG = pickle.load(open(f"{dataRoot}/PickleData/GraphsG.DHT.Data.p", "rb" ))
ethG = pickle.load(open(f"{dataRoot}/PickleData/GraphsG.EtOH.Data.p", "rb" ))


dhtBed = readBed(f"{dataRoot}/NodeBeds/creDHT.bed")
ethBed = readBed(f"{dataRoot}/NodeBeds/creEtOH.bed")



beds = {"dht": dhtBed, "eth": ethBed}

rangePlugIN(dhtG, beds)
rangePlugIN(ethG, beds)

print("Ranges are added!")





dhtFasta = readFasta(f"{dataRoot}/Sequences/cre.DHT.Seq.fasta")
ethFasta = readFasta(f"{dataRoot}/Sequences/cre.EtOH.Seq.fasta")

fastas = {"dht": dhtFasta, "eth": ethFasta}

seqPlugIN(dhtG, fastas)
seqPlugIN(ethG, fastas)

print("Sequences are added!")


GroDF = pickle.load(open(f"{dataRoot}/PickleData/Gro.DF.p", "rb" ))
for node in ethG.nodes():
    ethG.nodes[node]["lvlM"] = GroDF.loc[node, "dmso.-"]
    ethG.nodes[node]["lvlP"] = GroDF.loc[node, "dmso.+"]
    #
    ethG.nodes[node]["D"] = GroDF.loc[node, "Directionality"]

for node in dhtG.nodes():
    dhtG.nodes[node]["lvlM"] = GroDF.loc[node, "dht.-"]
    dhtG.nodes[node]["lvlP"] = GroDF.loc[node, "dht.+"]
    #
    dhtG.nodes[node]["D"] = GroDF.loc[node, "Directionality"]


pickle.dump(ethG,open(f"{dataRoot}/PickleData/GraphsG.EtOH.Data.p", "wb" ))

pickle.dump(dhtG,open(f"{dataRoot}/PickleData/GraphsG.DHT.Data.p", "wb" ))


# arpBed = readBed(f"{dataRoot}/Regions/ProteinBounds/AR.bed")
#
# ariBed = readBed(f"{dataRoot}/Regions/ProteinBounds/ARID.bed")
#
# ctcBed = readBed(f"{dataRoot}/Regions/ProteinBounds/CTCF.bed")
#
# ezhBed = readBed(f"{dataRoot}/Regions/ProteinBounds/EZH2.bed")
#
# foxBed = readBed(f"{dataRoot}/Regions/ProteinBounds/FOXA1.bed")
#
# medBed = readBed(f"{dataRoot}/Regions/ProteinBounds/MED1.bed")
#
# nmyBed = readBed(f"{dataRoot}/Regions/ProteinBounds/NMYC.bed")
#
# polBed = readBed(f"{dataRoot}/Regions/ProteinBounds/POLR2.bed")


#
# ProteinBeds = {"arp":arpBed,
#                "ctc":ctcBed,
#                "fox":foxBed,
#                "med":medBed,
#                "ezh":ezhBed,
#                "ari":ariBed,
#                "pol":polBed,
#                "nmy":nmyBed
#
# }
#
#
# featurePlugIN(G, ProteinBeds, "BoundProteins")
#
# # featurePlugIN(G, HistoneBeds, "Histones")
#
#
#
#
# GCmetBed = {}
# with open(f"{dataRoot}/DNAmet/GCint.bed") as tsvfile:
#     reader = csv.reader(tsvfile, delimiter='\t')
#     for row in reader:
#         # row3 = row[3].split(".")[0]
#         row3 = row[3]
#         GCmetBed[row3] = (row[0], int(row[1]), int(row[2]), int(row[4]), int(row[5]))
#
# GCmetBed = sortBed(GCmetBed)
# GCmetBed = {"met": GCmetBed}
#
#
# featurePlugIN(G, GCmetBed, "DNAmet", DNAmet=True)
#
#
#
#
# pickle.dump(C,open(f"{dataRoot}/PickleData/GraphsCData.p", "wb" ))
# pickle.dump(T,open(f"{dataRoot}/PickleData/GraphsTData.p", "wb" ))
pickle.dump(G,open(f"{dataRoot}/PickleData/GraphsGData.p", "wb" ))
