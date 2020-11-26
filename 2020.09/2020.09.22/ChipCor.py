


from Functions.Packages import *
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

chipAll = pickle.load(open(f"{dataRoot}/ChipAll.p", "rb" ))

chipAll2 = chipAll.set_index("Unnamed: 0").drop(columns=["nodeClass"])

#
# splitChip = chipAll.set_index("Unnamed: 0").drop(columns=["nodeClass"]).to_dict('split')
#
# # list(itertools.product(splitChip["index"], splitChip["index"])))
#
# con = chipAll[chipAll["nodeClass"] == "con"].set_index("Unnamed: 0").drop(columns=["nodeClass"]).to_dict('split')
# ind = chipAll[chipAll["nodeClass"] == "ind"].set_index("Unnamed: 0").drop(columns=["nodeClass"]).to_dict('split')
#


# len(list(itertools.product(con["index"], ind["index"])))


corrMat = chipAll2.T.corr()
pickle.dump(corrMat, open(f"{dataRoot}/CorrMatrixChip.p", "wb"))

mat = chipAll2.values.T

samples = chipAll2.shape[0]
centered_x = mat - np.sum(mat, axis=0, keepdims=True) / samples
centered_y = mat - np.sum(mat, axis=0, keepdims=True) / samples
cov_xy = 1./(samples - 1) * np.dot(centered_x.T, centered_y)
var_x = 1./(samples - 1) * np.sum(centered_x**2, axis=0)
var_y = 1./(samples - 1) * np.sum(centered_y**2, axis=0)
corrcoef_xy = cov_xy / np.sqrt(var_x[:, None] * var_y[None,:])



lut = dict(zip(chipAll["nodeClass"].unique(), ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]))
row_colors = chipAll["nodeClass"].map(lut)
DF_corr = pd.DataFrame(corrcoef_xy)

DF_dism = 1 - DF_corr   # distance matrix
linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')


fig = plt.figure(figsize=[10, 10])

# sns.clustermap(DF, row_colors=row_colors, row_linkage=linkage)
sns.clustermap(DF_corr, row_colors=row_colors, col_colors=row_colors, col_cluster=False)
# sns.heatmap(DF_corr)
fig.savefig(f"{figureRoot}/CorrMat3.pdf")
