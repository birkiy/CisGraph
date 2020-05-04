
from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))



intraCHRset = set()
for t in C.nodes():
    s = C.nodes()[t]["subG"]
    for e in list(s.edges()):
        intraCHRset.add(e)
        intraCHRset.add((e[1], e[0]))
Gset = set()
for e in list(G.edges()):
    Gset.add(e)
    Gset.add((e[1], e[0]))
interCHRset = Gset.difference(intraCHRset)



intraCHRset2 = set()
for e in intraCHRset:
    if not (e[1], e[0]) in intraCHRset2:
        intraCHRset2.add(e)
interCHRset2 = set()
for e in interCHRset:
    if not (e[1], e[0]) in interCHRset2:
        interCHRset2.add(e)
Gset2 = set()
for e in Gset:
    if not (e[1], e[0]) in Gset2:
        Gset2.add(e)





interCHR = {key: [] for key in list(itertools.combinations_with_replacement(["hom", "upP", "dwP", "con", "ind", "non"], 2))}
intraCHR = {key: [] for key in list(itertools.combinations_with_replacement(["hom", "upP", "dwP", "con", "ind", "non"], 2))}

for intra in intraCHRset2:
    edgeType = list(G.edges()[intra]["edgeType"])[0]
    edgeType2 = (edgeType[1], edgeType[0])


    if not edgeType in intraCHR.keys():
        intraCHR[edgeType2] += [intra]
    else:
        intraCHR[edgeType] += [intra]

for inter in interCHRset2:
    edgeType = list(G.edges()[inter]["edgeType"])[0]
    edgeType2 = (edgeType[1], edgeType[0])

    if not edgeType in interCHR.keys():
        interCHR[edgeType2] += [inter]
    else:
        interCHR[edgeType] += [inter]


print("interCHR", [len(_) for _ in interCHR.values()])
print("intraCHR", [len(_) for _ in intraCHR.values()])



edTy = list(itertools.combinations_with_replacement(["upP", "dwP", "con", "ind", "non"], 2))


interC = []
intraC = []

for ty in edTy:

    if "upP" in ty or "dwP" in ty:
        continue

    if (ty in interCHR.keys()) and (ty in intraCHR.keys()):
        interC += [len(set(interCHR[ty]))]
        intraC += [len(set(intraCHR[ty]))]





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

        if not ty in interCHR.keys():
            if len(set(interCHR[ty2])) == 0 or len(set(intraCHR[ty2])) == 0  :
                ratio = 0
                sum = 0
            else:
                ratio =( len(set(interCHR[ty2])) / len(set(intraCHR[ty2])) )
                sum = len(set(interCHR[ty2])) + len(set(intraCHR[ty2]))
        else:
            if len(set(interCHR[ty])) == 0 or len(set(intraCHR[ty])) == 0:
                ratio = 0
                sum = 0
            else:
                ratio =( len(set(interCHR[ty])) / len(set(intraCHR[ty])) )
                sum = len(set(interCHR[ty])) + len(set(intraCHR[ty]))
        R[i] += [ratio]
        S[i] += [sum]

    i += 1



upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "upP"]
dwP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "dwP"]

print(len(upP), len(dwP))


print(pd.DataFrame(R, nodeClasses, nodeClasses))
print(pd.DataFrame(S, nodeClasses, nodeClasses))
