


from Functions.Helpers import *


def fromGI(G, file, colorPalette):
    """
    fromGI() function takes a converted (please check the format) GI object and initialize a networkx graph.

        - Requires an empty nx graph object, converted file location and a color colorPalette (see formats)
        - Returns number of each node class.

    """

    nodeClassC = {}
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            abNode = row[0]

            aNode = abNode.split(",")[0]
            aNodeClass = aNode.split(".")[0]
            aNode = aNode.split(".",maxsplit=1)[1]

            bNode = abNode.split(",")[1]
            bNodeClass = bNode.split(".")[0]
            bNode = bNode.split(".",maxsplit=1)[1]

            weight = row[1]
            fdr = row[2]

            weight = int(weight)
            fdr = float(fdr)

            G.add_node(aNode,
                       color=colorPalette[aNodeClass],
                       nodeClass=aNodeClass,
                       logFC=0,
                       Qval=0,
                       nodeRange=(),
                       index=0,
                       nodeName=aNode,
                       tad="",
                       com="",
                       chr="",
                       subG=None,
                       subT=None,
                       subM=None)
            G.add_node(bNode,
                       color=colorPalette[bNodeClass],
                       nodeClass=bNodeClass,
                       logFC=0,
                       Qval=0,
                       nodeRange=(),
                       index=0,
                       nodeName=bNode,
                       tad="",
                       com="",
                       chr="",
                       subG=None,
                       subT=None,
                       subM=None)

            if (aNode, bNode) in G.edges:
                w2 = G.edges[(aNode, bNode)]["weight"]
                weight += w2
                G.edges[(aNode, bNode)]["weight"] = weight

            else:
                G.add_edge(aNode, bNode,
                weight=weight,
                fdr=fdr,
                edgeType={(aNodeClass, bNodeClass) : 0},
                color="#888888",
                distance="")


            if (not aNodeClass in nodeClassC.keys()):
                nodeClassC[aNodeClass] = 1
            elif (not bNodeClass in nodeClassC.keys()):
                nodeClassC[bNodeClass] = 1
            elif aNodeClass in nodeClassC.keys():
                nodeClassC[aNodeClass] += 1
            elif bNodeClass in nodeClassC.keys():
                nodeClassC[bNodeClass] += 1

    return nodeClassC



def fromGIup(L, U, upLvl="tad",lwLvl="subG"):
    """
    This function connects differnt level of graphs with given combination upper level and lower level.

        - Levels are G -> T -> M -> C, based on genomic organization.
        - Each level represented by a graph object, and their nodes are DNA regions in genome.

        - While going up in the system corresponded subGraphs are entegrated to one upper level.
        - At the same time current lower level nodes get an attribution that their positions in upper levels.

        - For fully connected system, please run this code as followings:
            [0]: fromGIup(G, T, "tad", "subG")
            [1]: fromGIup(G, M, "com", "subG")
            [2]: fromGIup(G, C, "chr", "subG")

            [3]: fromGIup(T, M, "com", "subT")
            [4]: fromGIup(T, C, "chr", "subT")

            [5]: fromGIup(M, C, "chr", "subM")
        - Note that firstly each upper level in the system should get the lowest subGraph attribution, then system can be builded.

    """

    Lbed = dict(L.nodes(data="nodeRange"))

    Lbed = {k: v for k,v in sorted(Lbed.items(), key=lambda kv: kv[1][1])}

    Ubed = dict(U.nodes(data="nodeRange"))

    Ubed = {k: v for k,v in sorted(Ubed.items(), key=lambda kv: kv[1][1])}

    Unodes = list(U.nodes())

    for node in Unodes:

        nodeRange = Ubed[node]
        nodeClass = U.nodes[node]["nodeClass"]

        nUp = rangesFromUpperRange(nodeRange[0], nodeRange[1], nodeRange[2], [Lbed],
                               ["L"])

        n = [_[0] for _ in nUp]

        if len(n) == 0:
            U.remove_node(node)
            continue


        for nodeL in n:
            L.nodes[nodeL][upLvl] = node

        U.nodes[node][upLvl] = node

        l = L.subgraph(n).copy()

        U.nodes[node][lwLvl] = l
