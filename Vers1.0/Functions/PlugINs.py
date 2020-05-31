

from Functions.Packages import *
from Functions.Helpers import *



def colorEdgesPlugIN(G, colorPalette):
    """
    This function plugs edge color depends on the degree of the current edge

        - Normally edges are #888888 as default.
        - After this function they attributed according to higher degree node anchor.

    G: Networkx graph that have been run fromGI() function before.
    """

    for edge in G.edges():
        aNode = edge[0]
        bNode = edge[1]

        aNodeClass = G.nodes[aNode]["nodeClass"]
        bNodeClass = G.nodes[bNode]["nodeClass"]

        if G.degree[aNode] > G.degree[bNode]:
            G[aNode][bNode]["color"] = colorPalette[aNodeClass]
        else:
            G[aNode][bNode]["color"] = colorPalette[bNodeClass]




def rangePlugIN(G, beds):
    """
    This function plugs the genomic positions of DNA element to the graph

        - Nodes will get "nodeRange" attribution of 3 element tuple (chr, start, end)
        - Edges will get "distance" attribution according to two anchor kb distance on genome, a string of a float.
        - Distance will get "INF" if they are at different chromosomes.

    G: Networkx graph that have been run fromGI() function before.
    beds: A python dictionary that contains nodeClasses as keys() and corresponding bed dictionaries as values(), created by readBed() function
    """

    for node in G.nodes():
        nodeClass = G.nodes[node]["nodeClass"]
        nodeRange = beds[nodeClass][node]

        G.nodes[node]["nodeRange"] = nodeRange


    for edge in G.edges():
        aNode = edge[0]
        bNode = edge[1]

        aNodeRange = G.nodes[aNode]["nodeRange"]
        bNodeRange = G.nodes[bNode]["nodeRange"]

        if aNodeRange[0] == bNodeRange[0]:
            distance = aNodeRange[1] - bNodeRange[1]
            distance = str(distance)
        else:
            distance = "INF"

        G.edges[(aNode, bNode)]["distance"] = distance


def enhancerFunctionPlugIN(G, Beds, colorPalette):

    for enh in [_ for _ in G.nodes() if G.nodes[_]["nodeClass"] == "enh"]:
        rangeX = G.nodes[enh]["nodeRange"]

        beds = list(Beds.values())
        nodeClasses = list(Beds.keys())
        nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2],beds ,nodeClasses )

        if len(nTuple) > 0:
            for node in nTuple:
                G.nodes[enh]["nodeClass"] = node[1]
                G.nodes[enh]["nodeName"] = node[0]
                G.nodes[enh]["color"] = colorPalette[G.nodes[enh]["nodeClass"]]



def logFCPlugIN(G, file, thQ=0.05, thLFC=(-1, +1), geneClass="pro", colorPalette={"upP": "#F95355", "dwP": "#4C6472"}):
    """
    This function plugs the logFC and Qval attribution to gene nodes in order to assign additional gene node classes "upP" and "dwP".

    G: Networkx graph that have been run fromGI() function before.
    file: CSV  format DEG analysis result file location,
        - Index[1]: geneSymbol, Index[3]: logFC, Index[7]: Qval
    thQ: Threshold for q value, default is 0.05
    thLFC: Threshold for log Fold-Change, default is (-1, +1)
    geneClass: First assigned node class for genes in G, default is "hom"
    colorPalette: A colorPalette (see formats) for newly assigned gene classes, default is {"upP": "#F95355", "dwP": "#4C6472"}
    """
    geneNodes = [node for node in G.nodes() if G.nodes[node]["nodeClass"] == geneClass]

    with open(file) as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        next(reader)
        for row in reader:
            geneSymbol = row[1]
            logFC = row[3]
            Qval = row[7]

            if not geneSymbol in geneNodes:
                continue

            if (logFC == "NA" or Qval == "NA"):
                continue

            logFC = float(logFC)
            Qval = float(Qval)



            if (logFC > thLFC[1] and Qval < thQ):
                nodeClass = list(colorPalette.keys())[0]
                G.nodes[geneSymbol]["color"] = colorPalette[nodeClass]
                G.nodes[geneSymbol]["nodeClass"] = nodeClass
                edges = list(G.edges(geneSymbol))
                for edge in edges:
                    aClass = G.nodes[edge[0]]["nodeClass"]
                    bClass = G.nodes[edge[1]]["nodeClass"]
                    G.edges[edge]["edgeType"] = {(aClass, bClass): 0}

            elif (logFC < thLFC[0] and Qval < thQ):
                nodeClass = list(colorPalette.keys())[1]
                G.nodes[geneSymbol]["color"] = colorPalette[nodeClass]
                G.nodes[geneSymbol]["nodeClass"] = nodeClass
                edges = list(G.edges(geneSymbol))
                for edge in edges:
                    aClass = G.nodes[edge[0]]["nodeClass"]
                    bClass = G.nodes[edge[1]]["nodeClass"]
                    G.edges[edge]["edgeType"] = {(aClass, bClass): 0}







def powerNodePlugIN(G, P, Beds, colorPalette):

            Gbed = dict(G.nodes(data="nodeRange"))

            Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}

            for pNode in Gbed:
                beds = list(Beds.values())
                nodeClasses = list(Beds.keys())
                rangeX = Gbed[pNode]
                nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2], beds, nodeClasses)
                if len(nTuple) < 1:
                    G.remove_node(pNode)
                    continue

                pn = nx.Graph()
                for node in nTuple:
                    if not node in P.nodes():

                        P.add_node(node[0],
                                    nodeClass=node[1],
                                    color=colorPalette[node[1]],
                                    elm=[pNode])

                        pn.add_node(node[0],
                                    nodeClass=node[1],
                                    color=colorPalette[node[1]],
                                    elm=[pNode])

                    else:
                        P.nodes[node]["elm"] += [pNode]

                        pn = G.nodes[pNode]["subP"].copy()

                        pn.nodes[node]["elm"] += [pNode]



                G.nodes[pNode]["subP"] = pn
