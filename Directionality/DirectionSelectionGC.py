



home = "/home/birkiy/dataGC/"



conFasta = readFasta(home + "Fasta/ARBS/cons-enc.fasta")
indFasta = readFasta(home + "Fasta/ARBS/ind.enc.fasta")
nonFasta = readFasta(home + "Fasta/ARBS/nonactive-enc.fasta")


# Note to store bed files as data frane

conBed = readBed(home + "cons-arbs.bed")
conBed = sortBed(conBed)
#
# conBed = pd.read_csv(home + "cons-arbs.bed", sep="\t", names=["#chrom", "start", "end", "name"])
# conBed["NodeClass"] = "con"
#
# indBed = pd.read_csv(home + "ind-arbs.bed", sep="\t", names=["#chrom", "start", "end", "name"])
# indBed["NodeClass"] = "ind"
#
# nonBed = pd.read_csv(home + "Non-Active-ARBS.bed", sep="\t", names=["#chrom", "start", "end", "name"])
# nonBed["NodeClass"] = "non"
#
indBed = readBed(home + "ind-arbs.bed")
indBed = sortBed(indBed)

nonBed = readBed(home + "Non-Active-ARBS.bed")
nonBed = sortBed(nonBed)

ArbsFasta = {"non":nonFasta,
             "ind": indFasta,
             "con": conFasta}

ArbsBed = {"non":nonBed,
             "ind": indBed,
             "con": conBed}

nodeClass = ("non", "ind", "con")
for node in nodeClass:
    arbsFasta = ArbsFasta[node]
    arbsBed = ArbsBed[node]
    arbsMbed = {}
    arbsPbed = {}
    for arb in arbsFasta.keys():
        seqStr = arbsFasta[arb]
        centerIdx = int(len(seqStr)/2)
        firstPart, secondPart = seqStr[:centerIdx], seqStr[centerIdx:]
        if GC(firstPart) > GC(secondPart):
            #Minus direction: End - 350
            rangeX = arbsBed[arb]
            arbsMbed[arb] = (rangeX[0], rangeX[1], rangeX[2])
        if GC(firstPart) < GC(secondPart):
            #Plus direction: Start + 350
            rangeX = arbsBed[arb]
            arbsPbed[arb] = (rangeX[0], rangeX[1], rangeX[2])
    writeBed(arbsMbed, node + ".+.bed")
    writeBed(arbsPbed, node + ".-.bed")

OvohF5el
writeBed(arbsMbed, "ARBS.+.bed")
writeBed(arbsPbed, "ARBS.-.bed")




arbsBeds = {"con":conBed,
               "ind":indBed,
               "non":nonBed
}


arbsFasta = {"con":conFasta,
               "ind":indFasta,
               "non":nonFasta
}


for node in arbsBeds.keys():
    arbsMbed = {}
    arbsPbed = {}
    for arb in arbsFasta[node]:
        seqStr = arbsFasta[node][arb]
        centerIdx = int(len(seqStr)/2)
        firstPart, secondPart = seqStr[:centerIdx], seqStr[centerIdx:]
        if GC(firstPart) >= GC(secondPart):
            #Minus direction: End - 350
            rangeX = arbsBeds[node][arb]
            arbsMbed[arb] = (rangeX[0], rangeX[1], rangeX[2]-350)
        if GC(firstPart) <= GC(secondPart):
            #Plus direction: Start + 350
            rangeX = arbsBeds[node][arb]
            arbsPbed[arb] = (rangeX[0], rangeX[1]+350, rangeX[2])
    # OvohF5el
    writeBed(arbsMbed, node + ".ARBS.M.bed")
    writeBed(arbsPbed, node + ".ARBS.P.bed")
