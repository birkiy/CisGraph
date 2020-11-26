


from Functions.Packages import *

from liftover import get_lifter

GR1i = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="1hr_induced")
GR1i["Condition"] = "1"
GR1i["Situation"] = "Induced"
GR1r = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="1hr_repressed")
GR1r["Condition"] = "1"
GR1r["Situation"] = "Repressed"

GR4i = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="4hr_induced")
GR4i["Condition"] = "4"
GR4i["Situation"] = "Induced"
GR4r = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="4hr_repressed")
GR4r["Condition"] = "4"
GR4r["Situation"] = "Repressed"

GR8i = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="8hr_induced")
GR8i["Condition"] = "8"
GR8i["Situation"] = "Induced"
GR8r = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="8hr_repressed")
GR8r["Condition"] = "8"
GR8r["Situation"] = "Repressed"

GR12i = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="12hr_induced")
GR12i["Condition"] = "12"
GR12i["Situation"] = "Induced"
GR12r = pd.read_excel("/home/ualtintas/GR/41467_2018_7607_MOESM3_ESM.xlsx", sheet_name="12hr_repressed")
GR12r["Condition"] = "12"
GR12r["Situation"] = "Repressed"


GR = pd.concat([GR1i, GR1r, GR4i, GR4r, GR8i, GR8r, GR12i, GR12r])


new = GR["coordinate"].str.split(":", expand=True)

GR["chr"] = new[0]

new2 = new[1].str.split("-", expand=True)
GR["start"] = new2[0]
GR["end"] = new2[1]

GR = GR.sort_values("start").reset_index(drop=True)



converter = get_lifter('hg38', 'hg19')
chrom = 'chr1'
pos = 103786442
converter[chrom][pos]
converter.query(chrom, pos)

GR['startHg19'] = GR.apply(lambda x: str(converter.query(x['chr'],int(x['start']))[0][1]) if len(converter.query(x['chr'],int(x['start']))) > 0 else None ,axis=1)
GR['endHg19'] = GR.apply(lambda x: str(converter.query(x['chr'],int(x['end']))[0][1]) if len(converter.query(x['chr'],int(x['end']))) > 0 else None ,axis=1)

AR = pd.read_csv("/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/DEG/lncap-merged-results-deseq2.txt", sep="\t")


# CPM stands for enhancer signal from starrseq and FC stands for fold change with DEX treatment
# ind => higher enhancer signal and higher foldchange with DEX treatment
# ind = GR[(GR["logCPM"] >= 1) & (GR["logFC"] >= 1)]
ind = GR[(GR["logFC"] >= 1)]
len(ind["coordinate"].unique()) # => 628
ind["nodeClass"] = "ind"


GR = GR[~(GR["logFC"] >= 1)]

# con => higher enhancer signal and lower foldchange, no change with DEX treatment
con = GR[(GR["logCPM"] >= 1) & ((GR["logFC"] < 1) & (GR["logFC"] > -1))]
len(con["coordinate"].unique()) # => 2070
con["nodeClass"] = "con"

GR = GR[~((GR["logCPM"] >= 1) & ((GR["logFC"] < 1) & (GR["logFC"] > -1)))]


# npn => lower enhancer signal and lower foldchange, no change with DEX treatment
non = GR[((GR["logCPM"] < 1) & (GR["logCPM"] > -1)) &  ((GR["logFC"] < 1) & (GR["logFC"] > -1))]
len(non["coordinate"].unique()) # => 5695
non["nodeClass"] = "non"

GR = pd.concat([ind, con, non])




fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

condition = 1

params = dict(
    x="logCPM",
    y="logFC",
    hue="nodeClass"
)


conditions = [1, 4, 8, 12]
for i, condition in enumerate(conditions):
    ax1 = fig.add_subplot(gs[i])
    sns.scatterplot(**params, data=GR[GR["Condition"] == str(condition)])
    plt.title(f"DEX treatment for {condition} hours")

fig.savefig(f"{figureRoot}/GRscatt.pdf")


# Condition problem
# coordinate     logFC    logCPM   nodeClass   condition
# chr10:100488703-100489413  1.040433  1.352457   ind 8
# chr10:100488703-100489413  0.399212  1.352457   con 4
# chr10:100488703-100489413  0.441258  1.352457   con 1
# chr10:100488703-100489413  0.759203  1.352457   con 12



sizes = [628, 2070, 5695]
labels = ["ind", "con", "non"]


fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(wspace=0.4)

conditions = [1, 4, 8, 12]
for i, condition in enumerate(conditions):
    sizes = [
        len(ind[ind["Condition"] == str(condition)]["coordinate"].unique()),
        len(con[con["Condition"] == str(condition)]["coordinate"].unique()),
        len(non[non["Condition"] == str(condition)]["coordinate"].unique())
        ]
    ax1 = fig.add_subplot(gs[i])
    plt.pie(sizes, labels=labels, autopct='%1.0f%%')
    plt.title(f"DEX treatment for {condition} hours\n(n={sum(sizes)})")

fig.savefig(f"{figureRoot}/GRenhPie.pdf")

#########################################3
GRr = GR[((GR["startHg19"] != None) | (GR["endHg19"] != None))]

bedGR = GRr.groupby(["chr", "startHg19", "endHg19","nodeClass"]).size().reset_index().drop(columns=0)

GR.groupby(["coordinate", "nodeClass"])

bedGR["name"] = [f"GRBS.Starr.{i}" for i in list(bedGR.index)]
conGR = bedGR[bedGR["nodeClass"] == "con"][["chr", "startHg19", "endHg19", "name"]]
conGR.to_csv('/home/ualtintas/GR/GR.con.bed', sep="\t", index=False, header=False)

nonGR = bedGR[bedGR["nodeClass"] == "non"][["chr", "startHg19", "endHg19", "name"]]
nonGR.to_csv('/home/ualtintas/GR/GR.non.bed', sep="\t", index=False, header=False)

indGR = bedGR[bedGR["nodeClass"] == "ind"][["chr", "startHg19", "endHg19", "name"]]
indGR.to_csv('/home/ualtintas/GR/GR.ind.bed', sep="\t", index=False, header=False)
