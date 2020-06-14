

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# BiocManager::install("CAGEr")


library(CAGEr)

inputFiles = list.files("/home/ualtintas/LncapExo/CAGE",
                        "processed.bam$",
                        full.names = TRUE)


inputFiles = c(inputFiles,
  c("/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.bam",
  "/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.bam")
)

library(BSgenome.Hsapiens.UCSC.hg19)


ce <- CAGEexp(genomeName     = "BSgenome.Hsapiens.UCSC.hg19" ,
              inputFiles     = inputFiles,
              inputFilesType = "bam" ,
              sampleLabels   = basename(inputFiles)
)


getCTSS(ce,
  removeFirstG=FALSE,
  sequencingQualityThreshold=-1,
  mappingQualityThreshold=30
)


# Plot the correlations
pdf(file = "/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage/cageCorPlot.pdf",
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches


corr.m <- plotCorrelation2( ce, samples = "all"
                          , tagCountThreshold = 1, applyThresholdBoth = FALSE
                          , method = "spearman")

dev.off()


# Plot the cumulatices
pdf(file = "/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage/cageCumPlot.pdf",
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)

dev.off()


mergeSamples(ce, mergeIndex = c(1,2,3,4,5,6,7,1,2,3,4,5,6,7),
            mergedSampleLabels = c("LNCaP.0h.rep1", "LNCaP.3h.rep1", "LNCaP.6h.rep1", "LNCaP.12h.rep1", "LNCaP.18h.rep1", "LNCaP.48h.rep1", "LNCaP.72h.rep1"))

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 150), alpha = 2.55, T = 1*10^5)




clusterCTSS( object = ce
         , threshold = 0.1
         , thresholdIsTpm = TRUE
         , nrPassThreshold = 1
         , method = "distclu"
         , maxDist = 20
         , removeSingletons = TRUE
         , keepSingletonsAbove = 5
)



samples = c("LNCaP.0h.rep1", "LNCaP.3h.rep1", "LNCaP.6h.rep1", "LNCaP.12h.rep1", "LNCaP.24h.rep1", "LNCaP.48h.rep1", "LNCaP.72h.rep1")

dht = tagClusters(ce, sample=samples[7])
eth = tagClusters(ce, sample=samples[1])

dhtP = dht[dht$strand %in% "+",]
dhtM = dht[dht$strand %in% "-",]





dhtTCPlu = dhtTC[dhtTC$strand %in% "+",]
dhtTCMin = dhtTC[dhtTC$strand %in% "-",]




exportCTSStoBedGraph(ce, values = "normalized", format = "BigWig")
exportCTSStoBedGraph(ce, values = "raw", format = "BigWig")


dhtTCMin[dhtTCMin$chr %in% "chr1",]
dhtTCPlu[dhtTCPlu$chr %in% "chr1",]



ethTC = tagClusters(ce, sample="cage.etoh.bam")




exportToBed(object = ce, what = "tagClusters",  oneFile = FALSE)


OvohF5el
sequencingQualityThreshold = -1


CAGEr::import.bam.ctss("/home/ualtintas/LncapExo/CAGE/cage.dht.bam", "bam" ,sequencingQualityThreshold = 10,
mappingQualityThreshold = 20)



grDHT = CAGEr:::import.bam.ctss("/home/ualtintas/LncapExo/CAGE/cage.dht.bam", filetype="bam" , sequencingQualityThreshold=-1, mappingQualityThreshold=30, removeFirstG=FALSE, genome="BSgenome.Hsapiens.UCSC.hg19")
grETH = CAGEr:::import.bam.ctss("/home/ualtintas/LncapExo/CAGE/cage.etoh.bam", filetype="bam" , sequencingQualityThreshold=-1, mappingQualityThreshold=30, removeFirstG=FALSE, genome="BSgenome.Hsapiens.UCSC.hg19")


CAGEr::bam2CTSS(grDHT, )



bam2CTSS
