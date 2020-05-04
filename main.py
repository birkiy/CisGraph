

print("\nGI code is running for G VCaP!\n")
subprocess.call (["/usr/bin/Rscript", "--vanilla", "GI2NX/GI.VCaP/GI.G.VCaP.R"])

print("\nGI code is running for T VCaP!\n")
subprocess.call (["/usr/bin/Rscript", "--vanilla", "GI2NX/GI.VCaP/GI.T.VCaP.R"])

print("\nGI code is running for C VCaP!\n")
subprocess.call (["/usr/bin/Rscript", "--vanilla", "GI2NX/GI.VCaP/GI.C.VCaP.R"])



print("\nGI Levels are initializing!\n")
import GI2NX.NetworkLevels.LevelInitGIs

print("\nGI Plugs!\n")
import GI2NX.NetworkLevels.LevelPlugINs

print("\nGI Levels are connecting!\n")
import GI2NX.NetworkLevels.LevelConnections
