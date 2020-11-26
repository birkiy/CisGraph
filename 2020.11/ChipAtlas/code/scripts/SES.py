
import multiprocessing
import numpy as np

# 2 bam files, IP, input
import deeptools
from deeptools import bamHandler
import deeptools.countReadsPerBin as countR

nThreads = 4

files = [
    "results/mapping/LNCaP.H3K4me3.etoh.4h.rep1.final.bam",
    "results/mapping/SRX4108929.1.control.final.bam"
]

# file = "results/mapping/SRX4108929.1.control.final.bam"
mappedReads = []
for file in files:
    mappedReads.append(bamHandler.openBam(file, returnStats=True, nThreads=nThreads)[1])


sizeFactorBasedOnMappedReads = np.array(mappedReads, dtype='float64')

sizeFactorBasedOnMappedReads = sizeFactorBasedOnMappedReads.min() / sizeFactorBasedOnMappedReads


cr = countR.CountReadsPerBin(files,
                                binLength=50,
                                numberOfSamples=10000,
                                extendReads=False,
                                numberOfProcessors=nThreads)





chromsizes, non_common = deeptools.utilities.getCommonChrNames([bamHandler.openBam(file) for file in files])


chrNames, chrLengths = list(zip(*chromsizes))

genomeSize = sum(chrLengths)


bam = bamHandler.openBam(file)


for read in bam.fetch("chr3", 9998999, 9999999):
    print(read.is_unmapped)


b = cr.run().transpose()



def calculateSESestimates(readCountsPerbin):
    if len(readCountsPerbin) != 2:
        readCountsPerbin = readCountsPerbin.transpose()
    p = np.sort(readCountsPerbin[0, :]).cumsum() # Chip
    q = np.sort(readCountsPerbin[1, :]).cumsum() # input
    # take the difference between max values and normalize accordingly
    diff = np.abs(p / p[-1] - q / q[-1])
    maxIndex = int(np.argmax(diff) * 0.8)
    cumSum = np.array(
        [
            float(p[maxIndex]),
            float(q[maxIndex])
            ]
        )
    sizeFactorsSES = cumSum.min() / cumSum
    meanSES = [
        np.mean(
            np.sort(
                readCountsPerbin[0, :]
                )[:maxIndex]
            ),
        np.mean(
            np.sort(
                readCountsPerbin[1, :]
                )[:maxIndex]
            )
        ]
    return sizeFactorsSES, meanSES



values = np.zeros(4)
maxNumReads = (np.percentile(readCountsPerbin[maxIndex:,0], 90))
values[0] = mean









while(maxIndex < len(p)):
    # in rare cases the maxIndex maps to a zero value.
    # In such cases, the next index is used until
    # a non zero value appears.
    print(maxIndex)
    cumSum = np.array([float(p[maxIndex]), float(q[maxIndex])])
    if cumSum.min() > 0:
        break
    maxIndex += 1





median = np.median(num_reads_per_bin, axis=1)

mean = []
std = []
for values in b:
    maxNumReads = (np.percentile(values, 90))
    if maxNumReads == 0:
        maxNumReads = (np.percentile(values, 99))
        if maxNumReads == 0:
            print("all genomic regions sampled from one ")
            "of the bam files have no reads.\n"
            values = values[values <= maxNumReads]
    mean.append(np.mean(values))
    std.append(np.std(values))


mean = np.array(mean)
readsPerBin = mean if avg_method == 'mean' else median
