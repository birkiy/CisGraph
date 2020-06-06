

from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))
ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))



colorPalette = ["#FADE89", "#57A4B1", "#B0D894"]

# Connected
b = nx.betweenness_centrality(G)
d = nx.degree_centrality(G)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])


B["Degree Centrality"] = D["Degree Centrality"]

# B["Degree Centrality"] += 1
# B["Betweenness Centrality"] += 1

nodeClasses = []
for node in B.index:
    nodeClasses += [G.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])

# Arbs["log(Betweenness Centrality)"] = np.log(Arbs["Betweenness Centrality"] + 1)

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")]




fig = plt.figure(figsize=[11,6])
gs = gridspec.GridSpec(ncols=2, nrows=1)

plt.subplots_adjust(hspace=0.1, wspace=0.2)

ax1 = fig.add_subplot(gs[0])

ax1 = sns.boxplot(x="nodeClass", y="Degree Centrality", data=Arbs, palette=colorPalette)
ax1.spines['bottom'].set_visible(False)
plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_ylim(0.003,0.01)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.title("Not Connected - Degree Centrality")
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
add_stat_annotation(ax1,x="nodeClass", y="Degree Centrality", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')

plt.title("Not Connected \n Degree Centrality")

ax2 = fig.add_subplot(gs[1])

ax = sns.boxplot(x="nodeClass", y="Degree Centrality", data=Arbs, palette=colorPalette, ax=ax2)
ax2.spines['top'].set_visible(False)
ax2.set_ylim(-10**(-8), 0.003)
ax2.set_ylabel("")
ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.text(0.08, 0.5, "Degree Centrality", va='center', rotation='vertical')



####

ax1 = fig.add_subplot(gs[1])

ax1 = sns.boxplot(x="nodeClass", y="Betweenness Centrality", data=Arbs, palette=colorPalette, ax=ax1, showfilers=False)
ax1.spines['bottom'].set_visible(False)
plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_ylim(1.02*10**(-5),0.03)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
add_stat_annotation(ax1,x="nodeClass", y="Betweenness Centrality", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')


plt.title("Not Connected \n Betweenness Centrality")


# ax2 = fig.add_subplot(gs[1,1])
# ax = sns.boxplot(x="nodeClass", y="Betweenness Centrality", data=Arbs, palette=colorPalette, ax=ax2)
# ax2.spines['top'].set_visible(False)
# ax2.set_ylim(-10**(-8), 1.02*10**(-5))
# ax2.set_ylabel("")
# ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#
#


B["Betweenness Centrality"][len(B["Betweenness Centrality"])/2]

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.text(0.5, 0.5, "Betweenness Centrality", va='center', rotation='vertical')


fig.savefig("Centrality.NotConnected.Seperated.pdf")


ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(dhtG))
}
gIdx = 0
tmpOldComp = []
for node in dhtG.nodes():
    tmpNewComp = set(nx.node_connected_component(dhtG, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]

# Largest => index 4
g = G.subgraph(ComponentsG[4])



b = nx.betweenness_centrality(g)
d = nx.degree_centrality(g)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])


B["Degree Centrality"] = D["Degree Centrality"]

nodeClasses = []
for node in B.index:
    nodeClasses += [g.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])



fig = plt.figure(figsize=[11,6])
gs = gridspec.GridSpec(ncols=2, nrows=2)

plt.subplots_adjust(hspace=0.1, wspace=0.2)

ax1 = fig.add_subplot(gs[0,0])

ax1 = sns.boxplot(x="nodeClass", y="Degree Centrality", data=Arbs, palette=colorPalette)
ax1.spines['bottom'].set_visible(False)
plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_ylim(0.015,0.04)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
add_stat_annotation(ax1,x="nodeClass", y="Degree Centrality", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')

plt.title("Connected \n Degree Centrality")

ax2 = fig.add_subplot(gs[1,0])

ax = sns.boxplot(x="nodeClass", y="Degree Centrality", data=Arbs, palette=colorPalette, ax=ax2)
ax2.spines['top'].set_visible(False)
ax2.set_ylim(0, 0.015)
ax2.set_ylabel("")
ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.text(0.07, 0.5, "Degree Centrality", va='center', rotation='vertical')



####

ax1 = fig.add_subplot(gs[0,1])

ax1 = sns.boxplot(x="nodeClass", y="Betweenness Centrality", data=Arbs, palette=colorPalette, ax=ax1, showbox=False, showcaps=False, showmeans=False, linewidth=0)
ax1.spines['bottom'].set_visible(False)
plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_ylim(0.075,0.5)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
add_stat_annotation(ax1,x="nodeClass", y="Betweenness Centrality", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')


plt.title("Connected \n Betweenness Centrality")



ax2 = fig.add_subplot(gs[1,1])
ax = sns.boxplot(x="nodeClass", y="Betweenness Centrality", data=Arbs, palette=colorPalette, ax=ax2)
ax2.spines['top'].set_visible(False)
ax2.set_ylim(0, 0.075)
ax2.set_ylabel("")
ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))





d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

fig.text(0.5, 0.5, "Betweenness Centrality", va='center', rotation='vertical')






fig.savefig("Centrality.Connected.ylim.pdf")






colorPalette = ["#FADE89", "#57A4B1", "#B0D894"]


colorPalette = (("non", "#B0D894"),
                ("ind", "#57A4B1"),
                ("con", "#FADE89"))


fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

Stat = ("Degree Centrality", "Betweenness Centrality")


for node, color in colorPalette:
    sns.distplot(Arbs[Arbs["nodeClass"] == node]["Betweenness Centrality"], color=color, hist=False)

colorPalette = ["#FADE89", "#57A4B1", "#B0D894"]
sns.pairplot(Arbs, hue="nodeClass", palette=colorPalette)

ax = fig.add_subplot(gs[0])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Degree Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    ax.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend(loc="lower right")
ax.set_ylabel("CDF")
plt.title("Not Connected \n Degree Centrality")
ax.set_xlabel("Degree Centrality")
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))

ax = fig.add_subplot(gs[1])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Betweenness Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    ax.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend(loc="lower right")
ax.set_ylabel("CDF")
plt.title("Not Connected \n Betweenness Centrality")
ax.set_xlabel("Betweenness Centrality")
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
axins = ax.inset_axes([0.5, 0.3, 0.47, 0.47])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Betweenness Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    print(F2)
    axins.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

axins.set_xlim([ax.get_xlim()[0], 0.005])
axins.set_ylim([0.75, ax.get_ylim()[1]])
axins.set_xticks([])
axins.set_yticks([])
ax.indicate_inset_zoom(axins)



fig.savefig("ZoomedCDF.pdf")



##################



ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]

# Largest => index 4
g = G.subgraph(ComponentsG[4])



b = nx.betweenness_centrality(g)
d = nx.degree_centrality(g)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])


B["Degree Centrality"] = D["Degree Centrality"]

nodeClasses = []
for node in B.index:
    nodeClasses += [g.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])





fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

Stat = ("Degree Centrality", "Betweenness Centrality")




ax = fig.add_subplot(gs[0])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Degree Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    ax.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend(loc="lower right")
ax.set_ylabel("CDF")
plt.title("Connected \n Degree Centrality")
ax.set_xlabel("Degree Centrality")
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))

ax = fig.add_subplot(gs[1])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Betweenness Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    ax.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

plt.legend(loc="lower right")
ax.set_ylabel("CDF")
plt.title("Connected \n Betweenness Centrality")
ax.set_xlabel("Betweenness Centrality")
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
axins = ax.inset_axes([0.5, 0.3, 0.47, 0.47])
for nodeClass, color in colorPalette:
    Z = Arbs[Arbs["nodeClass"] == nodeClass]["Betweenness Centrality"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    print(F2)
    axins.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)

axins.set_xlim([ax.get_xlim()[0], 0.005])
axins.set_ylim([0.75, ax.get_ylim()[1]])
axins.set_xticks([])
axins.set_yticks([])
ax.indicate_inset_zoom(axins)



fig.savefig("ZoomedCDF.pdf")













# Connected
b = nx.betweenness_centrality(G)
d = nx.degree_centrality(G)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]

nodeClasses = []
for node in B.index:
    nodeClasses += [G.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")]




#######################################


fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)



fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#FADE89", "#57A4B1", "#B0D894"]
ax=pt.half_violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, bw = .2, cut = 0., scale = "area", width = .6, inner = None)
ax=sns.stripplot( x = dx, y = dy, data = Arbs, palette = colorPalette, edgecolor = "white", size = 3, jitter = 1, zorder = 0)
ax.legend(labels=["con", "ind", "non"], loc="upper left")
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = True, boxprops = {'facecolor':'none', "zorder":10}, showfliers=True, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)

plt.title("Not Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')

fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)

ax=pt.half_violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, bw = .2, cut = 0., scale = "area", width = .6, inner = None)
ax=sns.stripplot( x = dx, y = dy, data = Arbs, palette = colorPalette, edgecolor = "white", size = 3, jitter = 1, zorder = 0)
ax.legend(labels=["con", "ind", "non"], loc="upper left")
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = True, boxprops = {'facecolor':'none', "zorder":10}, showfliers=True, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)

plt.title("Not Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')

fig.savefig("NotConnected.RainCloud.pdf")


#
# colorPalette = ["#FADE89", "#57A4B1", "#B0D894"]
# dx = "nodeClass"; dy = "Betweenness Centrality"; ort = "h"; pal = sns.color_palette(n_colors=1); sigma = .2; df = Arbs
# ax=pt.RainCloud(x = dx, y = dy, data = df, bw = sigma, width_viol = .6, ax = ax, orient = ort, move = .2)
#


###################################3
ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]

# Largest => index 4
g = G.subgraph(ComponentsG[4])

b = nx.betweenness_centrality(g)
d = nx.degree_centrality(g)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]
nodeClasses = []
for node in B.index:
    nodeClasses += [g.nodes()[node]["nodeClass"]]

B["nodeClass"] = nodeClasses
Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])
###########################################


fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)



fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)

ax=pt.half_violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, bw = .2, cut = 0., scale = "area", width = .6, inner = None)
ax=sns.stripplot( x = dx, y = dy, data = Arbs, palette = colorPalette, edgecolor = "white", size = 3, jitter = 1, zorder = 0)
ax.legend(labels=["con", "ind", "non"], loc="upper left")
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = True, boxprops = {'facecolor':'none', "zorder":10}, showfliers=True, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)

plt.title("Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')


fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)

ax=pt.half_violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, bw = .2, cut = 0., scale = "area", width = .6, inner = None)
ax=sns.stripplot( x = dx, y = dy, data = Arbs, palette = colorPalette, edgecolor = "white", size = 3, jitter = 1, zorder = 0)
ax.legend(labels=["con", "ind", "non"], loc="upper left")
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = True, boxprops = {'facecolor':'none', "zorder":10}, showfliers=True, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)

plt.title("Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside')

fig.savefig("Connected.RainCloud.pdf")

TAM 1 Malım








bw = .2, cut = 0., , width = .6









# Connected
b = nx.betweenness_centrality(G)
d = nx.degree_centrality(G)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]

nodeClasses = []
for node in B.index:
    nodeClasses += [G.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")]



fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Not Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)



fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Not Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.savefig("NotConnected.Centrality.Violin.pdf")


###################################3
ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]

# Largest => index 4
g = G.subgraph(ComponentsG[4])

b = nx.betweenness_centrality(g)
d = nx.degree_centrality(g)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]
nodeClasses = []
for node in B.index:
    nodeClasses += [g.nodes()[node]["nodeClass"]]

B["nodeClass"] = nodeClasses
Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])
###########################################


fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)



fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)


fig.savefig("Connected.Centrality.Violin.pdf")
