
from Functions.Packages import *

def fromGI(G, file, nodeClasses, colorPalette):
    conBed = {}
    readBed(conBed, "Data/regions/cons-arbs.bed")

    nonBed = {}
    readBed(nonBed, "Data/regions/Non-Active-ARBS.bed")

    indBed = {}
    readBed(indBed, "Data/regions/ind-arbs.bed")

    homBed = {}
    readBed(homBed, "Data/DEG/promoters_ann_5kb.bed")

    nodeClassC = {}
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            aNode = row[0]
            aNodeClass = row[1]
            aLgFC = row[2]
            aQval = row[3]
            bNode = row[4]
            bNodeClass = row[5]
            bLgFC = row[6]
            bQval = row[7]
            weight = row[8]
            fdr = row[9]


            if (not aNodeClass in nodeClassC.keys()):
                nodeClassC[aNodeClass] = 1
            elif (not bNodeClass in nodeClassC.keys()):
                nodeClassC[bNodeClass] = 1
            elif aNodeClass in nodeClassC.keys():
                nodeClassC[aNodeClass] += 1
            elif bNodeClass in nodeClassC.keys():
                nodeClassC[bNodeClass] += 1


            aLgFC = float(aLgFC)
            bLgFC = float(bLgFC)

            aQval = float(aQval)
            bQval = float(bQval)

            weight = int(weight)


            if (aNodeClass == "hom" and aLgFC > 1 and aQval < 0.05):
                aNodeClass = "upP"
            elif (aNodeClass == "hom" and aLgFC < -1 and aQval < 0.05):
                aNodeClass = "dwP"
            if (bNodeClass == "hom" and bLgFC > 1  and bQval < 0.05):
                bNodeClass = "upP"
            elif (bNodeClass == "hom" and bLgFC < -1 and bQval < 0.05):
                bNodeClass = "dwP"


            if (aNodeClass == "hom" or bNodeClass == "hom"):
                continue



            if aNodeClass == "non":
                aNodeRange = nonBed[aNode.split(".")[1]]
            elif aNodeClass == "ind":
                aNodeRange = indBed[aNode.split(".")[1]]
            elif aNodeClass == "con":
                aNodeRange = conBed[aNode.split(".")[1]]
            elif aNodeClass == "upP" or aNodeClass == "dwP" or aNodeClass == "hom" :
                aNodeRange = homBed[aNode.split(".")[1]]

            if bNodeClass == "non":
                bNodeRange = nonBed[bNode.split(".")[1]]
            elif bNodeClass == "ind":
                bNodeRange = indBed[bNode.split(".")[1]]
            elif bNodeClass == "con":
                bNodeRange = conBed[bNode.split(".")[1]]
            elif bNodeClass == "upP" or bNodeClass == "dwP" or bNodeClass == "hom" :
                bNodeRange = homBed[bNode.split(".")[1]]

            if aNodeRange[0] == bNodeRange[0]:
                distance = aNodeRange[1] - bNodeRange[1]
                distance = str(distance)
            else:
                distance = "INF"

            G.add_node(aNode,
                       color=colorPalette[nodeClasses.index(aNodeClass)],
                       nodeClass=aNodeClass,
                       lgFC=aLgFC,
                       Qval=aQval,
                       nodeRange=aNodeRange,
                       index=0,
                       nodeName=aNode,
                       tad="",
                       com="",
                       chr="")
            G.add_node(bNode,
                       color=colorPalette[nodeClasses.index(bNodeClass)],
                       nodeClass=bNodeClass,
                       lgFC=bLgFC,
                       Qval=bQval,
                       nodeRange=bNodeRange,
                       index=0,
                       nodeName=bNode,
                       tad="",
                       com="",
                       chr="")



            if (aNode, bNode) in G.edges:

                w2 = G.edges[(aNode, bNode)]["weight"]
                weight += w2
                G.edges[(aNode, bNode)]["weight"] = weight

            else:
                G.add_edge(aNode, bNode,
                           weight=weight,
                           edgeType=(aNodeClass, bNodeClass),
                           color="#888888",
                           distance=distance
                           )



    for edge in G.edges():
        aNode = edge[0]
        bNode = edge[1]

        aNodeClass = G.nodes[aNode]["nodeClass"]
        bNodeClass = G.nodes[bNode]["nodeClass"]

        if G.degree[aNode] > G.degree[bNode]:
            G[aNode][bNode]["color"] = colorPalette[nodeClasses.index(aNodeClass)]
        else:
            G[aNode][bNode]["color"] = colorPalette[nodeClasses.index(bNodeClass)]



def readBed(bed, file):
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]), int(row[2]))






def fromGIup(L, U, file, upLvl, upBed, **kwargs):

    # L lower level graph,
    # U upper level graph

    Lbed = dict(L.nodes(data="nodeRange"))

    Lbed = {k: v for k,v in sorted(Lbed.items(), key=lambda kv: kv[1][1])}

    """
    if kwargs["L2"]:
        L2bed = dict(kwargs["l2"].nodes(data="nodeRange"))

        L2bed = {k: v for k, v in sorted(L2bed.items(), key=lambda kv: kv[1][1])}
    """

    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            abNode = row[0]
            aNode = abNode.split(",")[0]
            bNode = abNode.split(",")[1]

            weight = row[1]
            fdr = row[2]

            aNodeRange = upBed[aNode.split(".", maxsplit=1)[1]]
            bNodeRange = upBed[bNode.split(".", maxsplit=1)[1]]


            aTAD = rangesFromTadRange(aNodeRange[0], aNodeRange[1], aNodeRange[2], [Lbed],
                                   ["G"])
            bTAD = rangesFromTadRange(bNodeRange[0], bNodeRange[1], bNodeRange[2], [Lbed],
                                   ["G"])

            """
            if kwargs["L2"]:
                aChr = rangesFromTadRange(aNodeRange[0], aNodeRange[1], aNodeRange[2], [L2bed],
                                          ["G"])
                bChr = rangesFromTadRange(bNodeRange[0], bNodeRange[1], bNodeRange[2], [L2bed],
                                          ["G"])

                xa = [_[0] for _ in aChr]
                xb = [_[0] for _ in bChr]

                for nodeA, nodeB in zip(xa, xb):
                    L.nodes[nodeA][kwargs["u2lvl"]] = aNode
                    L.nodes[nodeB][kwargs["u2lvl"]] = bNode

                xga = L.subgraph(ca).copy()
                xgb = L.subgraph(cb).copy()
            """



            ca = [_[0] for _ in aTAD]
            cb = [_[0] for _ in bTAD]

            if len(ca) == 0 or len(cb) == 0:
                continue

            for nodeA, nodeB in zip(ca, cb):
                L.nodes[nodeA][upLvl] = aNode
                L.nodes[nodeB][upLvl] = bNode

            ga = L.subgraph(ca).copy()
            gb = L.subgraph(cb).copy()

            # Add subG sub T

            U.add_node(aNode,
                       color="#888888",
                       nodeClass=upLvl,
                       Qval=fdr,
                       nodeRange=aNodeRange,
                       index=0,
                       nodeName=aNode,
                       subG=ga)
            U.add_node(bNode,
                       color="#888888",
                       nodeClass=upLvl,
                       Qval=fdr,
                       nodeRange=bNodeRange,
                       index=0,
                       nodeName=bNode,
                       subG=gb)

            aNodeClass = U.nodes[aNode]["nodeClass"]
            bNodeClass = U.nodes[bNode]["nodeClass"]

            U.add_edge(aNode, bNode,
                       weight=weight,
                       edgeType={(aNodeClass, bNodeClass) : 0},
                       color="#888888")












def rangesFromTadRange(chrKey, start, end, beds, bedname):
    values = []
    for idx, enhancer in enumerate(beds):
        enhancerEnd = [(value[2], value[1], value[0], key) for key, value in enhancer.items()]
        left = bi.bisect_left(enhancerEnd, (start, end))

        enhancerStart = [(value[1], value[2], value[0], key) for key, value in enhancer.items()]
        right = bi.bisect_right(enhancerStart, (end, start))
        values += [
            (keys, bedname[idx])
            for keys, value in list(enhancer.items())[left:right]
            if chrKey == value[0] and
               value[2] - start > start - value[1] and
               end - value[1] > value[2] - end
        ]
    return values
