

from Functions.GIFunctions import *

C = pickle.load(open(f"{dataRoot}/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(f"{dataRoot}/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


print("Level connection starts!")

print("G -> T")
fromGIup(G, T, "tad", "subG")
# print("G -> M")
# fromGIup(G, M, "com", "subG")
print("G -> C")
fromGIup(G, C, "chr", "subG")

# print("T -> M")
# fromGIup(T, M, "com", "subT")
print("T -> C")
fromGIup(T, C, "chr", "subT")

# print("M -> C")
# fromGIup(M, C, "chr", "subM")



print("Level connection is done! Updated numbes:")

print("You have an G level graph of %i nodes, %i edges, %i components." % (
    len(G.nodes),
    len(G.edges),
    nx.number_connected_components(G)
))

print("You have an T level graph of %i nodes, %i edges, %i components." % (
    len(T.nodes),
    len(T.edges),
    nx.number_connected_components(T)
))

print("You have an C level graph of %i nodes, %i edges, %i components." % (
    len(C.nodes),
    len(C.edges),
    nx.number_connected_components(C)
))



pickle.dump(C,open(f"{dataRoot}/tmpData/GraphsCData.p", "wb" ))
pickle.dump(T,open(f"{dataRoot}/tmpData/GraphsTData.p", "wb" ))
pickle.dump(G,open(f"{dataRoot}/tmpData/GraphsGData.p", "wb" ))
