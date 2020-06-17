

from Functions.Packages import *

# subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/Vers2.0/GI2NX/GI.VCaP/GI.VCaP.R"])

subprocess.call ([envR , "--vanilla", f"{projectRoot}/GI2NX/GI.HiC/GI.HiC.R"])

# import GI2NX.BuildNX.InitializeGraph
# print("Init done!\n")


#import GI2NX.BuildNX.RegulatoryElementPlugIns
#print("PlugIn done!\n")


# import GI2NX.BuildNX.ConnectGraphs
# print("Connect done!\n")
