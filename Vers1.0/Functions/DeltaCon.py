

from Functions.Packages import *

def generateD(A):
    n,m = A.shape
    degs = A.sum(axis=1).flatten()
    D = sps.spdiags(degs,[0],n,n)
    return D

def generateS(A, e=10**-5):
    D = generateD(A)
    I = sps.identity(A.shape[0])
    S = (I + (e**2)*D - e*A)**-1
    return S

def RootED(S1, S2):
    RED = 0
    n, m = S1.shape
    for i in range(n):
        for j in range(m):
            RED += (math.sqrt(S1[i,j]) - math.sqrt(S2[i,j]))**2
    d = math.sqrt(RED)
    return d
    # d = math.sqrt(((np.sqrt(S1[:,1])-np.sqrt(S2[:,1]))**2).sum())

def sim(d):
    return 1/(1+d)

def deltaCon0(A1, A2):
    S1 = generateS(A1)
    S2 = generateS(A2)
    return sim(RootED(S1,S2))

def deltaConAttrNode(A1, A2):
    # INPUT: affinity matrices S1, S2
    # edge files of G1(V, E1) and G2(V, E2), i.e., A1 and A2
    S1 = generateS(A1)
    S2 = generateS(A2)
    n = A1.shape[0]
    w = []
    for v in range(n):
        # If an edge adjacent to the node has changed, the node is responsible
        if sum(abs(A1[:,v] - A2[:,v])) > 0:
            w.append(RootED(S1[:,v], S2[:,v]))
        else:
            w.append(0)
    return w

    
def deltoConAtrrEdge(A1, A2, w):
    n = A1.shape[0]
    for v in range(n):
        srcW = w[v]
        r = A2[:,v] - A1[:,v]
        for k in range(n):
            destNode = k
            destW = w[k]
            if r[k] == 1:
                edgeScore = srcW + destW
                E.append(srcW, destNode, "+", edgeScore)
            if r[k] == -1:
                edgeScore = srcW + destW
                E.append(srcW, destNode, "-", edgeScore)
    return E


def NX2deltaCon(G1, G2, returnG=False):
    nodesU = set(G1.nodes()).union(set(G2.nodes()))

    for node in nodesU:
        if not node in G1.nodes():
            att = G2.nodes[node]
            G1.add_node(node)
            for at in att:
                G1.nodes[node][at] = att[at]
        if not node in G2.nodes():
            att = G1.nodes[node]
            G2.add_node(node)
            for at in att:
                G2.nodes[node][at] = att[at]

    _G1 = nx.Graph()
    _G2 = nx.Graph()
    for node in nodesU:
        _G1.add_node(node, nodeClass=G1.nodes[node]["nodeClass"], w= None)
        _G2.add_node(node, nodeClass=G2.nodes[node]["nodeClass"], w= None)

    _G1.add_edges_from(list(G1.edges()))
    _G2.add_edges_from(list(G2.edges()))

    A1 = nx.to_numpy_matrix(_G1)
    A2 = nx.to_numpy_matrix(_G2)

    d = deltaCon0(A1, A2)
    print(f"DeltaCon Similarity score of two graphs is {sim(d)}")

    w = deltaConAttrNode(A1, A2)
    print(f"Top 20 node impact {w}")
    a = sorted([(idx, wi) for idx, wi in zip(range(len(w)), w)], key=lambda kv: kv[1])
    a = np.array([a]*len(a))
    #
    A1 = np.array(list(map(lambda x, y: y[x], np.argsort(a), A1)))
    A2 = np.array(list(map(lambda x, y: y[x], np.argsort(a), A2)))
    print("Adjacency matrices are sorted!")
    E = deltaConAttrEdge(A1, A2, sorted(w))

    if returnG:
        for node, wi in zip(nodesU, w):
            _G1.nodes[node]["w"] = wi
            _G2.nodes[node]["w"] = wi
        return d, w, E, _G1, _G2

    return d, w, E
