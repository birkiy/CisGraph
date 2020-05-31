
from Functions.FeaturePlugIns import *


C = pickle.load(open(home + "/Data/PickleData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/PickleData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/PickleData/GraphsGData.p", "rb" ))



tadBed = readBed(home + "/Data/Regions/TAD.C42B.bed")

chrBed = readBed(home + "/Data/Regions/hg19.Chr.bed")

creBed = readBed(home + "/Data/Regions/creDHS.bed")


beds = {"cre": creBed, "tad": tadBed, "chr": chrBed}

rangePlugIN(G, beds)
rangePlugIN(T, beds)
rangePlugIN(C, beds)
print("Ranges are added!")





creFasta = readFasta(home + "/Data/Sequences/creDHS.fasta")

fastas = {"cre": creFasta}

seqPlugIN(G, fastas)

print("Sequences are added!")



arpBed = readBed(home + "/Data/Regions/ProteinBounds/AR.bed")

ariBed = readBed(home + "/Data/Regions/ProteinBounds/ARID.bed")

ctcBed = readBed(home + "/Data/Regions/ProteinBounds/CTCF.bed")

ezhBed = readBed(home + "/Data/Regions/ProteinBounds/EZH2.bed")

foxBed = readBed(home + "/Data/Regions/ProteinBounds/FOXA1.bed")

medBed = readBed(home + "/Data/Regions/ProteinBounds/MED1.bed")

nmyBed = readBed(home + "/Data/Regions/ProteinBounds/NMYC.bed")

polBed = readBed(home + "/Data/Regions/ProteinBounds/POLR2.bed")



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
with open(home + "/Data/DNAmet/GCint.bed") as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        # row3 = row[3].split(".")[0]
        row3 = row[3]
        GCmetBed[row3] = (row[0], int(row[1]), int(row[2]), int(row[4]), int(row[5]))

GCmetBed = sortBed(GCmetBed)
GCmetBed = {"met": GCmetBed}


featurePlugIN(G, GCmetBed, "DNAmet", DNAmet=True)




pickle.dump(C,open(home + "/Data/PickleData/GraphsCData.p", "wb" ))
pickle.dump(T,open(home + "/Data/PickleData/GraphsTData.p", "wb" ))
pickle.dump(G,open(home + "/Data/PickleData/GraphsGData.p", "wb" ))
