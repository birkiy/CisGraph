


from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


H = G.to_directed()

D = nx.DiGraph()
for edge in H.edges():
    d = G.edgeG = nx.DiGraph()s
