

from Functions.Helpers import *


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




def seqPlugIN(G, fastas):

    for node in G.nodes():
        nodeClass = G.nodes[node]["nodeClass"]
        nodeSequence = fastas[nodeClass][node]

        G.nodes[node]["nodeSequence"] = nodeSequence

        G.nodes[node]["nodeLength"] = len(nodeSequence)

        G.nodes[node]["GC"] = GC(nodeSequence)


def featurePlugIN(G, ProteinBeds, feature, DNAmet=False):
    Gbed = dict(G.nodes(data="nodeRange"))
    Gbed = sortBed(Gbed)
    if DNAmet:
        DNAmetBeds = ProteinBeds
        DNAmetMean, DNAmetLen = rangeMapping(Gbed, DNAmetBeds, DNAmet=True)
        DNAmetLen = {k:v for k,v in sorted(DNAmetLen.items(), key=lambda kv: kv[0])}
        DNAmetMean = {k:v for k,v in sorted(DNAmetMean.items(), key=lambda kv: kv[0])}
        for node in DNAmetMean.keys():
            G.nodes[node]["DNAmet"] = (DNAmetMean[node], DNAmetLen[node])
    else:
        Proteins, mappingsProteins = rangeMapping(Gbed, ProteinBeds)
        for node in mappingsProteins.keys():
            G.nodes[node][feature] = mappingsProteins[node]
