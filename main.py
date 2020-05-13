
from Functions.Packages import *


if conditionCheck["createGIs"] == "VCaP":
    print("\nGI code is running for G VCaP!\n")
    subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.G.VCaP.R"])

    print("\nGI code is running for T VCaP!\n")
    subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.T.VCaP.R"])

    print("\nGI code is running for C VCaP!\n")
    subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/GI.VCaP/GI.C.VCaP.R"])
elif conditionCheck["createGIs"] == "dhs":
    print("\nGI code is running for G VCaP - DHS!\n")
    subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/GI2NX/NDRmodel/GI.G.VCaP.R"])

else:
    print("GI creation is passed!")



if conditionCheck["changeGIs"] == "VCaP":
    print("Take VCaP")
    subprocess.call (["/bin/cp", home + "/Data/GIs/VCaP/*", "/Data/GIs/."])
elif conditionCheck["changeGIs"] == "HiChip":
    print("Take HiChip")
    subprocess.call (["/bin/cp", home + "/Data/GIs/HiChip/*", "/Data/GIs/."])
elif conditionCheck["changeGIs"] == "dhs":
    print("Take dhs")
    subprocess.call (["/bin/cp", home + "/Data/GIs/dhs/*", "/Data/GIs/."])


if conditionCheck["levelInit"]:
    print("\nGI Levels are initializing!\n")
    import GI2NX.NetworkLevels.LevelInitGIs
else:
    print("Level initiation is passed!")

if conditionCheck["levelPlugIN"]:
    print("\nGI Plugs!\n")
    import GI2NX.NetworkLevels.LevelPlugINs
else:
    print("Level plugINs is passed!")


if conditionCheck["levelConnect"]:
    print("\nGI Levels are connecting!\n")
    import GI2NX.NetworkLevels.LevelConnections
else:
    print("Level connection is passed!")



print("\nInter TADs!\n")
# import InterConnections.InterTADs

print("\nInter TADs Heatmap!\n")
# import InterConnections.Heatmap

print("\nY Chromosome!\n")
# import PlotGraphs.Ychrom

print("\nAll Chromosomes!\n")
# import PlotGraphs.Chr2Tad

# import PlotGraphs.TADneighbors

if conditionCheck["plotTadBFS"]:
    import PlotGraphs.TadBFSplots
    print("\nBFS plotted!\n")


if conditionCheck["motifSource"]:
    import MotifSearch.DaisyChains
    print("\nMotif plotted!\n")
