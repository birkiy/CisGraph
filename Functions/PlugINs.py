

from Functions.Packages import *



def colorEdgesPlugIN(G):
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
            G[aNode][bNode]["color"] = colorPalette[nodeClasses.index(aNodeClass)]
        else:
            G[aNode][bNode]["color"] = colorPalette[nodeClasses.index(bNodeClass)]




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




def logFCPlugIN(G, file, thQ=0.05, thLFC=(-1, +1), geneClass="hom", colorPalette={"upP": "#F95355", "dwP": "#4C6472"}):
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
        for row in reader:
            geneSymbol = row[1]
            logFC = row[3]
            Qval = row[7]

            logFC = float(logFC)
            Qval = float(Qval)

            if not geneSymbol in geneNodes:
                continue

            if (logFC > thLFC[1] and Qval < thQ):
                nodeClass = "upP"
                G.nodes[geneSymbol]["color"] = colorPalette["upP"]
            elif (logFC < thLFC[0] and Qval < thQ):
                nodeClass = "dwP"
                G.nodes[geneSymbol]["color"] = colorPalette["dwP"]
            else:
                nodeClass = geneClass


            for node in geneNodes:
                G.nodes[node]["nodeClass"] = nodeClass
                G.nodes[node]["logFC"] = logFC
                G.nodes[node]["Qval"] = Q
