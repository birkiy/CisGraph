
import os

home = "/home/birkiy/github/CisGraph/Vers2.0"

os.chdir(home)

env = "/home/birkiy/anaconda3/envs/CisGraph"

from Functions.Packages import *

# subprocess.call ([env + "/bin/Rscript", "--vanilla", "/home/birkiy/github/CisGraph/Vers2.0/GI2NX/GI.VCaP/GI.VCaP.R"])


# import GI2NX.BuildNX.InitializeGraph
# print("Init done!\n")


#import GI2NX.BuildNX.RegulatoryElementPlugIns
#print("PlugIn done!\n")


import GI2NX.BuildNX.ConnectGraphs
print("Connect done!\n")
