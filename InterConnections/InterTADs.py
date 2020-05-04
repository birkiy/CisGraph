
from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))



intraTADset = set()
for t in T.nodes():
    s = T.nodes()[t]["subG"]
    for e in list(s.edges()):
        intraTADset.add(e)
        intraTADset.add((e[1], e[0]))
Gset = set()
for e in list(G.edges()):
    Gset.add(e)
    Gset.add((e[1], e[0]))
interTADset = Gset.difference(intraTADset)



intraTADset2 = set()
for e in intraTADset:
    if not (e[1], e[0]) in intraTADset2:
        intraTADset2.add(e)
interTADset2 = set()
for e in interTADset:
    if not (e[1], e[0]) in interTADset2:
        interTADset2.add(e)
Gset2 = set()
for e in Gset:
    if not (e[1], e[0]) in Gset2:
        Gset2.add(e)





interTAD = {key: [] for key in list(itertools.combinations_with_replacement(["hom", "upP", "dwP", "con", "ind", "non"], 2))}
intraTAD = {key: [] for key in list(itertools.combinations_with_replacement(["hom", "upP", "dwP", "con", "ind", "non"], 2))}

for intra in intraTADset2:
    edgeType = list(G.edges()[intra]["edgeType"])[0]
    edgeType2 = (edgeType[1], edgeType[0])


    if not edgeType in intraTAD.keys():
        intraTAD[edgeType2] += [intra]
    else:
        intraTAD[edgeType] += [intra]

for inter in interTADset2:
    edgeType = list(G.edges()[inter]["edgeType"])[0]
    edgeType2 = (edgeType[1], edgeType[0])

    if not edgeType in interTAD.keys():
        interTAD[edgeType2] += [inter]
    else:
        interTAD[edgeType] += [inter]


print("interTAD", [len(_) for _ in interTAD.values()])
print("intraTAD", [len(_) for _ in intraTAD.values()])



edTy = list(itertools.combinations_with_replacement(["upP", "dwP", "con", "ind", "non"], 2))


interT = []
intraT = []

for ty in edTy:

    if "upP" in ty or "dwP" in ty:
        continue

    if (ty in interTAD.keys()) and (ty in intraTAD.keys()):
        interT += [len(set(interTAD[ty]))]
        intraT += [len(set(intraTAD[ty]))]





nodeClasses = ["hom", "upP", "dwP",
               "con", "ind", "non"]




i = 0
R = []
S = []
for row in nodeClasses:
    R += [[]]
    S += [[]]
    for col in nodeClasses:
        ty = (row, col)
        ty2 = (col, row)

        if not ty in interTAD.keys():
            if len(set(interTAD[ty2])) == 0 or len(set(intraTAD[ty2])) == 0  :
                ratio = 0
                sum = 0
            else:
                ratio =( len(set(interTAD[ty2])) / len(set(intraTAD[ty2])) )
                sum = len(set(interTAD[ty2])) + len(set(intraTAD[ty2]))
        else:
            if len(set(interTAD[ty])) == 0 or len(set(intraTAD[ty])) == 0:
                ratio = 0
                sum = 0
            else:
                ratio =( len(set(interTAD[ty])) / len(set(intraTAD[ty])) )
                sum = len(set(interTAD[ty])) + len(set(intraTAD[ty]))
        R[i] += [ratio]
        S[i] += [sum]

    i += 1



upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "dwP"]

print(len(upP), len(dwP))



print(pd.DataFrame(R, nodeClasses, nodeClasses))
print(pd.DataFrame(S, nodeClasses, nodeClasses))
