

from Functions.Helpers import *
from Functions.Packages import *



home = "/home/birkiy/github/CisGraph"



conBed = {}
readBed(conBed, home + "/Data/Regions/cons-arbs.bed")
conBed = {k: v for k,v in sorted(conBed.items(), key=lambda kv: kv[1][1])}

nonBed = {}
readBed(nonBed, home + "/Data/Regions/Non-Active-ARBS.bed")
nonBed = {k: v for k,v in sorted(nonBed.items(), key=lambda kv: kv[1][1])}

indBed = {}
readBed(indBed, home + "/Data/Regions/ind-arbs.bed")
indBed = {k: v for k,v in sorted(indBed.items(), key=lambda kv: kv[1][1])}






colorPalette ={"hom": "#888888",
               "con": "#FADE89",
               "ind": "#57A4B1",
               "non": "#B0D894",
               "tad": "#4c5172",
               "chr": "#5a4c72"
               }





homBed = {}
readBed(homBed, home + "/Data/Regions/promoters_ann_5kb.bed")

tadBed = {}
readBed(tadBed, home + "/Data/Regions/TAD.C42B.bed")

chrBed = {}
readBed(chrBed, home + "/Data/Regions/hg19.Chr.bed")

file = home + "/Data/development/EnhancerDomain.bed"

with open(file) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        chrK = row[0]
        start = int(row[1])
        end = int(row[2])

        enhID = row[3]


        nTuple = rangesFromUpperRange(chrK, start, end, [conBed, indBed, nonBed], ["con", "ind", "non"])

        pn = nx.Graph()
        for node in nTuple:
            pn.add_node(node[0],
                        nodeClass=node[1],
                        color=colorPalette[node[1]])

        G.nodes[enhID]["subP"] = pn
