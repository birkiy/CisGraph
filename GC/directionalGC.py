


# Script Location: CisGraph/GC/directionalGC.py
# Utils: Vers2.0

from GC.Utils import *

home = "/home/birkiy/dataGC/"


tsPfasta = readFasta(home + "Fasta/tsP.fasta")
tsMfasta = readFasta(home + "Fasta/tsM.fasta")

firstP = []
secondP = []
for dnaKey in tsPfasta:
    seqStr = tsPfasta[dnaKey]
    centerIdx = int(len(seqStr)/2)
    firstPart, secondPart = seqStr[:centerIdx], seqStr[centerIdx:]
    firstP += [GC(firstPart)]
    secondP += [GC(secondPart)]

firstM = []
secondM = []
for dnaKey in tsMfasta:
    seqStr = tsMfasta[dnaKey]
    centerIdx = int(len(seqStr)/2)
    firstPart, secondPart = seqStr[:centerIdx], seqStr[centerIdx:]
    firstM += [GC(firstPart)]
    secondM += [GC(secondPart)]





fig = plt.figure(figsize=(6,10))
gs = gridspec.GridSpec(ncols=2, nrows=2, height_ratios=[9,1])
# plt.subplots(sharey=True)

fig.add_subplot(gs[1,:])
gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
hc = ['#FF9154', '#FFB58F', '#FFFFFF', '#D6F5EC', '#58D0A6']
th = [0, 0.2, 0.5, 0.8, 1]
cdict = NonLinCdict(th, hc)
cm = LinearSegmentedColormap('test', cdict)

plt.imshow(gradient, aspect='auto', cmap=cm)
plt.axis("off")
plt.xlabel("Regulatory Element")

yDec = 0.04
y1 = 1.1
ylim = [0.2,1.2]

data = {"FirstPart": firstP,
        "SecondPart": secondP
        }
colorPalette = ["#FFA96B", "#75E3C1"]
ax1 = fig.add_subplot(gs[0,0])

boxPlot(data=data, combs=[[0,1]], colorPalette=colorPalette, y1=y1, yDec=yDec, ylim=ylim, paired=True)
plt.ylabel("GC")

data = {"FirstPart": firstM,
        "SecondPart": secondM
        }

ax2 = fig.add_subplot(gs[0,1], sharey=ax1)

boxPlot(data=data, combs=[[0,1]], colorPalette=colorPalette, y1=y1, yDec=yDec, ylim=ylim, paired=True)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax1, title='Plus Strand')
plt.setp(ax2, title='Minus Strand')
plt.setp(ax2, ylabel="")

fig.savefig("PromoterDirections.pdf")

# plt.show()
