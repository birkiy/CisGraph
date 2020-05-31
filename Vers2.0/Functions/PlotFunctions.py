
from Functions.Utils import *


def NonLinCdict(steps, hexcol_array):
    cdict = {'red': (), 'green': (), 'blue': ()}
    for s, hexcol in zip(steps, hexcol_array):
        rgb =matplotlib.colors.hex2color(hexcol)
        cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
        cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
        cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)
    return cdict



def boxPlot(data, combs, colorPalette, yDec, y1, ylim, multiTest=False, paired=False):
    PVal =[]
    for c in combs:
        if paired:
            pVal = scipy.stats.wilcoxon(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
        else:
            pVal = scipy.stats.ranksums(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
        PVal += [pVal]
    if multiTest:
        pAdj = multipletests(PVal, method='bonferroni')
    else:
        pAdj = PVal
    cGC = flat([data[d] for d in data])
    nodeClasses = [nodeClass for nodeClass in data.keys() for _ in range(len(data[nodeClass]))]
    df2 = {
        "nodeClass": nodeClasses,
        "GC": cGC
    }
    df2 = pd.DataFrame(df2)
    sns.boxplot(x="nodeClass", y="GC", data=df2, palette=colorPalette)
    # sns.swarmplot(x="nodeClass", y="GC", data=df2, color=".25")
    plt.xlabel("")
    plt.ylabel("")
    plt.ylim(ylim)
    for c, q in zip(combs, pAdj):
        pVal = q
        x1 = c[0]
        x2 = c[1]
        if y1 != 0.98:
            y1 = y1 - yDec
        y2 = y1 + (yDec / 3)
        if pVal < 0.0001:
            pS = "***"
        elif pVal < 0.001:
            pS = "**"
        elif pVal < 0.01:
            pS = "*"
        elif pVal < 0.05:
            pS = "."
        else:
            pS = "ns"
        plt.plot([x1, x1, x2, x2], [y1, y2, y2, y1], linewidth=1, color='k')
        plt.text((x1 + x2) * .5, y2, pS, ha='center', va='bottom', color="k", fontsize=10)
