


import pandas as pd
from snakemake.io import expand, glob_wildcards

sampleDF = pd.read_table("code/samples.tsv")

controlDF = pd.read_table("code/controls.tsv")

srrDF = pd.concat(
    [
        controlDF[["SRR", "Library"]],
        sampleDF[["SRR", "Library"]]
    ]
)

srxDF = pd.concat(
    [
        controlDF[["SRX", "SRR", "Library"]],
        sampleDF[["SRX", "SRR", "Library"]]
    ]
)





bigWigsMerged = [
    f"results/bigwig/genomcov/merged/{raw}.bigWig"
    for raw in set(sampleDF["SampleName"])
]


plotCorrelation = [
    f"results/plots/correlation.{samples}.{corr}.{plot}.pdf"
    for samples in ["final", "merged"]
    for corr in ["pearson", "spearman"]
    for plot in ["heatmap"]
]

# plotFingerprint = [
#     f"results/plots/fingerprint.{samples}.pdf"
#     for samples in ["rep", "merged"]
# ]

# desiredOutputList = bigWigs + plotCorrelation + plotFingerprint + bigWigsMerged
desiredOutputList = bigWigsMerged + plotCorrelation





def getLib(wildcards):
        """
        This function is used by MappingBowtie2's getFastqs and Bigwig rules' getParamsBamCom, getParamsBamCov helper functions.
        if wildcards has length 3 it came from bigwig rule otherwise mapping.
        """
        # print(wildcards, "getlib wildcards")
        # print(type(wildcards), len(wildcards))
        if wildcards.raw.find("control") > 0:
                # print("Enter getlib Control")
                raw = wildcards.raw.rsplit(".", 1)[0]
                if len(wildcards) == 3:
                        # print("Enter getlib Control stage")
                        if wildcards.stage == "merged":
                                # print("Enter getlib Control stage merged")
                                lib = controlDF.loc[controlDF["ControlName"] == raw, "Library"].to_list()[0]
                        else:
                                # print("Enter getlib Control stage merged else")
                                lib = controlDF.loc[controlDF["Raw"] == raw, "Library"].to_list()[0]
                else:
                        # print("Enter getlib Control stage else")
                        lib = controlDF.loc[controlDF["Raw"] == raw, "Library"].to_list()[0]

        else:
                # print("Enter getlib Samples")
                if len(wildcards) == 3:
                        # print("Enter getlib Samples stage")
                        if wildcards.stage == "merged":
                                # print("Enter getlib Samples stage merged")
                                lib = sampleDF.loc[sampleDF["SampleName"] == wildcards.raw, "Library"].to_list()[0]
                        else:
                                # print("Enter getlib Samples stage merged else")
                                lib = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "Library"].to_list()[0]
                else:
                        # print("Enter getlib Samples stage else")
                        lib = sampleDF.loc[sampleDF["Raw"] == wildcards.raw, "Library"].to_list()[0]
        return lib
