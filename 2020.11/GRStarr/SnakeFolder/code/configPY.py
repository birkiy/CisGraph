

import pandas as pd


sampleDF = pd.read_table("code/samples.tsv", dtype=str)
sampleDF = sampleDF.sort_index()


sampleDF["sampleName"] = sampleDF[["tf", "duration" , "treatment"]].agg(".".join, axis=1)

sampleDF.loc[sampleDF["duration"] == "input", "sampleName"] = "input.GSE114063"

#
# plots = [
#     f"coverage/plot.Heatmap.RPKM.pdf",
#     f"coverage/plot.Heatmap.BPM.pdf",
#     f"coverage/plot.Heatmap.None.pdf",
#     f"coverage/plot.Heatmap.SES.pdf",
#     f"coverage/plot.Heatmap.read.pdf"
# ]


uniqueFragment = [
    f"coverage/unique.library.pool{num}.bedpe"
    for num in range(1,13)
]

coverage = [
    f"coverage/countTableRaw.merged.txt",
    f"coverage/countTableRaw.rep.txt"
]

# desiredOutputList = coverage + ["mapping/processed/input.GSE114063.merged.final.bam"] #+ plots
desiredOutputList = uniqueFragment
