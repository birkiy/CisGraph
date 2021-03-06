

from Functions.Packages import *

def generateD(A):
    n,m = A.shape
    degs = A.sum(axis=1).flatten()
    D = sps.spdiags(degs,[0],n,n)
    return D

def generateS(A, e=None):
    D = generateD(A)
    if e is not None:
        e = e
    else:
        e = 1/(1+D.todense().max())
    I = sps.identity(A.shape[0])
    S = (I + (e**2)*D - e*A)**-1
    return S

def RootED(S1, S2, psa=False, c=None):
    RED = 0
    n, m = S1.shape
    if psa:
        if c == 0:
            d = math.sqrt(psa[n,c])
        else:
            d = math.sqrt(psa[n,c] - psa[n,c-1])
    else:
        for i in range(n):
            for j in range(m):
                RED += (math.sqrt(S1[i,j]) - math.sqrt(S2[i,j]))**2
        d = math.sqrt(RED)
    return d
    # d = math.sqrt(((np.sqrt(S1[:,1])-np.sqrt(S2[:,1]))**2).sum())


def RootEDpsa(S1, S2):
    n, m = S1.shape
    psa = np.array((n,m))
    for i in range(n):
        for j in range(m):
            a = ((math.sqrt(S1[i,j]) - math.sqrt(S2[i,j]))**2)
            if i == 0 and j == 0:
                psa[i][j] = a
            elif i == 0 and j > 0:
                psa[i][j] = psa[i][j-1] + a
            elif i > 0 and j == 0:
                psa[i][j] = psa[i-1][j] + a
            else:
                psa[i][j] = psa[i-1][j] + psa[i][j-1] - psa[i-1][j-1] + a

    return psa
    # d = math.sqrt(((np.sqrt(S1[:,1])-np.sqrt(S2[:,1]))**2).sum())







def sim(d):
    return 1/(1+d)

def deltaCon0(A1, A2, e=None):
    S1 = generateS(A1, e=e)
    S2 = generateS(A2, e=e)
    return sim(RootED(S1,S2))

def deltaConAttrNode(A1, A2, e=None):
    # INPUT: affinity matrices S1, S2
    # edge files of G1(V, E1) and G2(V, E2), i.e., A1 and A2
    S1 = generateS(A1, e=e)
    S2 = generateS(A2, e=e)
    n = A1.shape[0]
    w = []
    # psa = RootEDpsa(S1, S2)
    for v in range(n):
        # If an edge adjacent to the node has changed, the node is responsible
        if sum(abs(A1[:,v] - A2[:,v])) > 0:
            w.append(RootED(S1[:,v], S2[:,v]))
        else:
            w.append(0)
    return w


def deltaConAttrEdge(A1, A2, w):
    E = []
    n = A1.shape[0]
    for v in range(n):
        srcW = w[v]
        srcNode = v
        r = A2[:,v] - A1[:,v]
        for k in range(n):
            destNode = k
            destW = w[k]
            if r[k] == 1:
                edgeScore = srcW + destW
                E.append((srcNode, destNode, "+", edgeScore))
            if r[k] == -1:
                edgeScore = srcW + destW
                E.append((srcNode, destNode, "-", edgeScore))
    return E


def NX2deltaCon(G1, G2, returnG=False, e=None):
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


    d = deltaCon0(A1, A2, e=e)
    print(f"\nDeltaCon Similarity score of two graphs is {sim(d)}")

    w = deltaConAttrNode(A1, A2, e=e)
    print(f"\nTop 20 node impact {sorted(w, reverse=True)[:20]}")
    a = sorted([(idx, wi) for idx, wi in zip(range(len(w)), w)], key=lambda kv: kv[1])
    b = [_[0] for _ in a]
    a = np.array([b]*len(b))
    #
    A1 = np.array(list(map(lambda x, y: y[x], np.argsort(a), np.array(A1))))
    # np.array(list(map(lambda x, y: y[x], np.argsort(b), x)))
    A2 = np.array(list(map(lambda x, y: y[x], np.argsort(a), np.array(A2))))
    print("\nAdjacency matrices are sorted!")
    E = deltaConAttrEdge(A1, A2, sorted(w, reverse=True))

    if returnG:
        for node, wi in zip(nodesU, w):
            _G1.nodes[node]["w"] = wi
            _G2.nodes[node]["w"] = wi
        return d, w, E, _G1, _G2

    return d, w, E
