
import pandas as pd


sampleDF = pd.read_table("code/samples.tsv")

controlDF = pd.read_table("code/controls.tsv")

srrDF = pd.concat(
    [
        controlDF[["SRR", "Library"]],
        sampleDF[["SRR", "Library"]]
    ]
)



union = [
    f"results/coverage/{stage}/union/{fragment}.union"
    for stage in ["merged", "each"]
    for fragment in ["unique", "all"]
]


desiredOutputList = union

# Replicates
# uniques = [
#     f"results/coverage/unique.{raw}.final.bedpe"
#     for raw in set(list(sampleDF["Raw"]) + list(controlDF["Raw"] +".control"))
# ]
#
# bigWigs = [
#     f"results/bigwig/{type}/bamcov/final/{raw}.bigWig"
#     for raw in set(list(sampleDF["Raw"]))
#     for type in ["RPKM", "CPM", "unique"]
# ]
#
# # Samples
# uniquesMerged = [
#     f"results/coverage/unique.{raw}.merged.bedpe"
#     for raw in set(list(sampleDF["SampleName"]) + list(controlDF["ControlName"] +".control"))
# ]
# bigWigsMerged = [
#     f"results/bigwig/{type}/bamcom/merged/{raw}.bigWig"
#     for raw in set(sampleDF["SampleName"])
#     for type in ["RPKM", "CPM", "unique"]
# ]
#
#
# plotCorrelation = [
#     f"results/plots/correlation.{samples}.{corr}.{plot}.pdf"
#     for samples in ["final", "merged"]
#     for corr in ["pearson", "spearman"]
#     for plot in ["heatmap"]
# ]


# desiredOutputList = uniques + bigWigs + uniquesMerged + bigWigsMerged + plotCorrelation
#
# desiredOutputList = ["results/mapping/final/GSE114063.1.control.bam"
# #desiredOutputList = plotCorrelation




################################################################################

#
# plotFingerprint = [
#     f"results/plots/fingerprint.{samples}.pdf"
#     for samples in ["rep", "merged"]
# ]
#
#
#
# for i in range(100):
#     a = random.randint(1,10) # computer
#     b = random.randint(1,10) # mister
#     L = [random.randint(97,122) for _ in range(a)]
#     sL = [random.randint(97,122) for _ in range(b)]
#     alphabet = list(range(97,122))
#     tL = []
#     n = 0
#     while a > n:
#         check = alphabet[n]
#         if check in sL:
#             n += 1
#             a += 1
#             continue
#         else:
#             tL += [check]
#             n += 1
#     s = "".join([chr(l) for l in sL])
#     t = "".join([chr(l) for l in tL])
#     L += sL + tL
#     # print(f"s:{s} \t\t\t t:{t}")
#
#
# vals = input("Enter ints:")
#
#
#
# a, b, l, r = vals.split(" ")
#
#
#
# a = int(a)
# b = int(b)
# l = int(l)
# r = int(r)
# # a = random.randint(1,10) # computer
# # b = random.randint(1,10) # mister
# import time
#
# start = time.time()
# alphabet = list(range(97,122))
# a, b, l, r = [3,7,4,10**8]
# L = alphabet[:a]
# tL = alphabet[:a]
# while True:
#     a, b, l, r = [3,7,4,10**8]
#     sL = sorted([tL[-1] for _ in range(b)])
#     tL = []
#     n = 0
#     while a > n:
#         check = alphabet[n]
#         if check in sL:
#             n += 1
#             a += 1
#             continue
#         else:
#             tL += [check]
#             n += 1
#     s = "".join([chr(l) for l in sL])
#     t = "".join([chr(l) for l in tL])
#     L += sL + tL
#     if len(L) > l and len(L) > r:
#         S = set(L[l-1:r])
#         end = time.time()
#         print(len(S), start - end)
#         break
#
#
#
#
# print()
