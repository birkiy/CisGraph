



from Functions.Packages import *


logFC = pd.read_csv(f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv")
logFC = logFC[["Gene", "log2FoldChange", "padj"]]
logFC = logFC.sort_values("Gene").reset_index().drop("index", axis=1)
logFC["Gene"] = logFC["Gene"].astype(str)


logFC["-log(Qval)"] = logFC["padj"]


fig = plt.figure(figsize=(4, 5))
gs = gridspec.GridSpec(ncols=1, nrows=1)
plt.subplots_adjust(wspace=0.4)
