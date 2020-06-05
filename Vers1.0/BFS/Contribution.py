



for root in both:
    fig = plt.figure(figsize=[15,12])
    gs = gridspec.GridSpec(ncols=2, nrows=1)
    #
    axDht = fig.add_subplot(gs[0])
    gl = dhtRooted[root]
    sl = dhtShell[root]
    SingleGraphPlot(Gs[0], axDht, root, gl, sl)
    #
    axEth = fig.add_subplot(gs[1])
    gl = ethRooted[root]
    sl = ethShell[root]
    SingleGraphPlot(Gs[1], axEth, root, gl, sl)
    #
    fig.savefig(f"{figureRoot}/Cicero/{root}.pdf")
    plt.close("all")



Contributions = {}
for root in both:
    dl = dhtRooted[root]
    d = dhtG.subgraph(dl).copy()
    el = ethRooted[root]
    e = ethG.subgraph(el).copy()
    dU = set(dl).difference(set(el))
    contD = dict.fromkeys(dU, set())
    eU = set(el).difference(set(dl))
    contE = dict.fromkeys(eU, set())
    for u in dU:
        contD[u] = set(d.edges(u))
    for u in eU:
        contE[u] = set(e.edges(u))
    Contributions[root] = (contD, contE)


ConD = {"ind": {"nc": 0}, "con": {"nc": 0}, "non": {"nc": 0}}
for root in Contributions:
    nodes = Contributions[root][0]
    for n in nodes.keys():
        nc = dhtG.nodes[n]["nodeClass"]
        ConD[nc]["nc"] += 1
        edges = nodes[n]
        for ed in edges:
            edTy = list(dhtG.edges[ed]["edgeType"].keys())[0]
            if not edTy in ConD[nc].keys():
                ConD[nc][edTy] = 1
            ConD[nc][edTy] += 1


ConE = {"ind": {"nc": 0}, "con": {"nc": 0}, "non": {"nc": 0}}
for root in Contributions:
    nodes = Contributions[root][1]
    for n in nodes.keys():
        nc = ethG.nodes[n]["nodeClass"]
        ConE[nc]["nc"] += 1
        edges = nodes[n]
        for ed in edges:
            edTy = list(ethG.edges[ed]["edgeType"].keys())[0]
            if not edTy in ConE[nc].keys():
                ConE[nc][edTy] = 1
            ConE[nc][edTy] += 1
