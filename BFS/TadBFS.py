from Functions.CustomBFS import *



home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))

# Given list to tadL will create a dictionary that contains BFS tree node list and accorded shell structure with given depthLimit.
tadL = [G.nodes["STEAP4"]["tad"], G.nodes["KLK2"]["tad"], G.nodes["FKBP5"]["tad"]]

depthLimit = 2
tTree = {}

Shell = {}
RootedTAD = {}
for rootT in tadL:
    shell = [None for _ in range(depthLimit+1)]
    currentDepth = []
    root = rootT

    t = bfs_treeCustom(T, root, depth_limit=depthLimit, current_depth=currentDepth, distance=False,shell=shell)
    lvlN = depthLimit - min(currentDepth)
    if not lvlN in tTree:
        tTree[lvlN] = {}
        if not root in tTree[lvlN]:
            tTree[lvlN][root] = [t]
        else:
            tTree[lvlN][root] += [t]
    else:
        if not root in tTree[lvlN]:
            tTree[lvlN][root] = [t]
        else:
            tTree[lvlN][root] += [t]
    shell = [i for i in shell if i]
    # if root.find("KLK")+1 or root == "hom.FKBP5" or root == "hom.STEAP4" or root == "hom.ZBTB16":
    Shell[root] = shell
    RootedTAD[root] = list(t.nodes())
    #print(root, lvlN)
