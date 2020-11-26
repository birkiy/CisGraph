

library(ChIPpeakAnno)


home = "/home/ualtintas/github/Data/CisGraph/Vers1.0"

conFile = paste(home, "Regions/cons-arbs.bed", sep="/")
indFile = paste(home, "Regions/ind-arbs.bed", sep="/")
nonFile = paste(home, "Regions/Non-Active-ARBS.bed", sep="/")

con <- toGRanges(conFile, format="BED", header=FALSE)
ind <- toGRanges(indFile, format="BED", header=FALSE)
non <- toGRanges(nonFile, format="BED", header=FALSE)

ol <- findOverlapsOfPeaks(con, ind, non)

library(EnsDb.Hsapiens.v75)

conL = ol$peaklist["con"]
indL = ol$peaklist["ind"]
nonL = ol$peaklist["non"]

annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)


aCR<-assignChromosomeRegion(con, nucleotideLevel=TRUE,
                           precedence=c("Promoters", "immediateDownstream",
                                         "fiveUTRs", "threeUTRs",
                                         "Exons", "Introns"),
                           TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)

con.anno <- annotatePeak(con,
                                AnnotationData=annoData,
                                output="nearestBiDirectionalPromoters",
                                bindingRegion=c(-2000, 500))



annotatePeakInBatch(con, mart, featureType = c("TSS", "miRNA","Exon"),



pdf("/home/ualtintas/github/Figures/CisGraph/Vers1.0/ChipAnnot.pdf", width=8.5, height=5,onefile=T)

binOverFeature(conL, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count",
               main="Distribution of Constitutive peak numbers around TSS")


binOverFeature(indL, annotationData=annoData,
              radius=5000, nbins=20, FUN=length, errFun=0,
              ylab="count",
              main="Distribution of Inducible peak numbers around TSS")


binOverFeature(nonL, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count",
               main="Distribution of Inactive peak numbers around TSS")
