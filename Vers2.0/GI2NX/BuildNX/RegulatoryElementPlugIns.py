
from Functions.FeaturePlugIns import *


C = pickle.load(open(f"{dataRoot}/PickleData/GraphsCData.p", "rb" ))
T = pickle.load(open(f"{dataRoot}/PickleData/GraphsTData.p", "rb" ))
G = pickle.load(open(f"{dataRoot}/PickleData/GraphsGData.p", "rb" ))



tadBed = readBed(f"{dataRoot}/Regions/TAD.C42B.bed")

chrBed = readBed(f"{dataRoot}/Regions/hg19.Chr.bed")

creBed = readBed(f"{dataRoot}/Regions/creDHS.bed")


beds = {"cre": creBed, "tad": tadBed, "chr": chrBed}

rangePlugIN(G, beds)
rangePlugIN(T, beds)
rangePlugIN(C, beds)
print("Ranges are added!")





creFasta = readFasta(f"{dataRoot}/Sequences/creDHS.fasta")

fastas = {"cre": creFasta}

seqPlugIN(G, fastas)

print("Sequences are added!")



arpBed = readBed(f"{dataRoot}/Regions/ProteinBounds/AR.bed")

ariBed = readBed(f"{dataRoot}/Regions/ProteinBounds/ARID.bed")

ctcBed = readBed(f"{dataRoot}/Regions/ProteinBounds/CTCF.bed")

ezhBed = readBed(f"{dataRoot}/Regions/ProteinBounds/EZH2.bed")

foxBed = readBed(f"{dataRoot}/Regions/ProteinBounds/FOXA1.bed")

medBed = readBed(f"{dataRoot}/Regions/ProteinBounds/MED1.bed")

nmyBed = readBed(f"{dataRoot}/Regions/ProteinBounds/NMYC.bed")

polBed = readBed(f"{dataRoot}/Regions/ProteinBounds/POLR2.bed")



ProteinBeds = {"arp":arpBed,
               "ctc":ctcBed,
               "fox":foxBed,
               "med":medBed,
               "ezh":ezhBed,
               "ari":ariBed,
               "pol":polBed,
               "nmy":nmyBed

}


featurePlugIN(G, ProteinBeds, "BoundProteins")

# featurePlugIN(G, HistoneBeds, "Histones")




GCmetBed = {}
with open(f"{dataRoot}/DNAmet/GCint.bed") as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        # row3 = row[3].split(".")[0]
        row3 = row[3]
        GCmetBed[row3] = (row[0], int(row[1]), int(row[2]), int(row[4]), int(row[5]))

GCmetBed = sortBed(GCmetBed)
GCmetBed = {"met": GCmetBed}


featurePlugIN(G, GCmetBed, "DNAmet", DNAmet=True)




pickle.dump(C,open(f"{dataRoot}/PickleData/GraphsCData.p", "wb" ))
pickle.dump(T,open(f"{dataRoot}/PickleData/GraphsTData.p", "wb" ))
pickle.dump(G,open(f"{dataRoot}/PickleData/GraphsGData.p", "wb" ))
