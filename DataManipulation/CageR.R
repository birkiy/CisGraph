

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("CAGEr")


library(CAGEr)

inputFiles = c("/home/ualtintas/LncapExo/CAGE/cage.dht.bam",
               "/home/ualtintas/LncapExo/CAGE/cage.etoh.bam")


# inputFiles = list.files( system.file("extdata", package = "CAGEr")
#                        , "ctss$"
#                        , full.names = TRUE)

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


corr.m <- plotCorrelation2( ce, samples = "all"
                          , tagCountThreshold = 1, applyThresholdBoth = FALSE
                          , method = "spearman")



normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)

clusterCTSS( object = ce
         , threshold = 1
         , thresholdIsTpm = TRUE
         , nrPassThreshold = 1
         , method = "distclu"
         , maxDist = 20
         , removeSingletons = TRUE
         , keepSingletonsAbove = 5
)



dhtTC = tagClusters(ce, sample="cage.dht.bam")
dhtTCPlu = dhtTC[dhtTC$strand %in% "+",]
dhtTCMin = dhtTC[dhtTC$strand %in% "-",]




exportCTSStoBedGraph(ce, values = "normalized", format = "BigWig")



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
