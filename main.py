
from Functions.Packages import *

print("\nGI code is running for G VCaP!\n")
# subprocess.call (["/home/birkiy/anaconda3/envs/CisGraph/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.G.VCaP.R"])
print("Passed")

print("\nGI code is running for T VCaP!\n")
# subprocess.call (["/home/birkiy/anaconda3/envs/CisGraph/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.T.VCaP.R"])
print("Passed")

print("\nGI code is running for C VCaP!\n")
# subprocess.call (["/home/birkiy/anaconda3/envs/CisGraph/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.C.VCaP.R"])
print("Passed")



print("\nGI Levels are initializing!\n")
# import GI2NX.NetworkLevels.LevelInitGIs
print("Passed")

print("\nGI Plugs!\n")
# import GI2NX.NetworkLevels.LevelPlugINs
print("Passed")

print("\nGI Levels are connecting!\n")
# import GI2NX.NetworkLevels.LevelConnections
print("Passed")

print("\nInter TADs!\n")
# import InterConnections.InterTADs

print("\nInter TADs Heatmap!\n")
# import InterConnections.Heatmap

print("\nY Chromosome!\n")
# import PlotGraphs.Ychrom

print("\nAll Chromosomes!\n")
# import PlotGraphs.Chr2Tad

#import PlotGraphs.TADneighbors


import PlotGraphs.TadBFSplots
