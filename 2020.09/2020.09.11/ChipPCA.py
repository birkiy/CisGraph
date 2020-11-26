
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


from Functions.Packages import *

chipARBS = pd.read_csv(f"{dataRoot}/ARBS_chip_seq_2.csv")

chipProm = pd.read_csv(f"{dataRoot}/Promoters_chip_seq.csv")

chip = pd.concat([chipARBS, chipProm])



conBed = pd.read_csv(home + "/ARBSs/regions/cons-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
indBed = pd.read_csv(home + "/ARBSs/regions/ind-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
nonBed = pd.read_csv(home + "/ARBSs/regions/Non-Active-ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nARBed = pd.read_csv(home + "/ARBSs/regions/negativeControl.ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])

tsPBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.+.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])
tsMBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.-.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])

conChip = chip[chip["Unnamed: 0"].isin(conBed["name"])]
conChip["nodeClass"] = "con"
indChip = chip[chip["Unnamed: 0"].isin(indBed["name"])]
indChip["nodeClass"] = "ind"
nonChip = chip[chip["Unnamed: 0"].isin(nonBed["name"])]
nonChip["nodeClass"] = "non"
nARChip = chip[chip["Unnamed: 0"].isin(nARBed["name"])]
nARChip["nodeClass"] = "nAR"

tsPChip = chip[chip["Unnamed: 0"].isin(tsPBed["name"])]
tsPChip["nodeClass"] = "tss"
tsMChip = chip[chip["Unnamed: 0"].isin(tsMBed["name"])]
tsMChip["nodeClass"] = "tss"

chipAll = pd.concat([conChip, indChip, nonChip, nARChip, tsPChip, tsMChip])

chipAll = chipAll.reset_index(drop=True)



##################################################33


features = list(chipAll.columns)
features.remove("Unnamed: 0")
features.remove("nodeClass")

x = chipAll.loc[:, features].values

y = chipAll.loc[:, ["nodeClass"]].values

x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])


finalDf = pd.concat([principalDf, chipAll[['nodeClass']]], axis = 1)



fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = ['con', 'ind', 'non', "nAR", "tss"]
colors = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['nodeClass'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)


ax.legend(targets)
ax.grid()


fig.savefig(f"{figureRoot}/PCA.Chip.pdf")
