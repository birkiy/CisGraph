

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggupset)


.libPaths(c("/kuacc/users/ualtintas20/apps/R/lib64"))

home = "/kuacc/users/ualtintas20/github/Data/CisGraph/Vers1.0"
home = "/home/ualtintas/github/Data/CisGraph/Vers1.0"

conFile = paste(home, "Regions/cons-arbs.bed", sep="/")
indFile = paste(home, "Regions/ind-arbs.bed", sep="/")
nonFile = paste(home, "Regions/Non-Active-ARBS.bed", sep="/")
othFile = paste(home, "Regions/otherARBS.bed", sep="/")
tsPFile = "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed"
tsMFile = "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed"


conAnnFile = paste(home, "Regions/con.AnnFile.csv", sep="/")
indAnnFile = paste(home, "Regions/ind.AnnFile.csv", sep="/")
nonAnnFile = paste(home, "Regions/non.AnnFile.csv", sep="/")
othAnnFile = paste(home, "Regions/oth.AnnFile.csv", sep="/")
tsPAnnFile = "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.AnnFile.csv"
tsMAnnFile = "/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.AnnFile.csv"


con <- readPeakFile(conFile)
ind <- readPeakFile(indFile)
non <- readPeakFile(nonFile)
oth <- readPeakFile(othFile)
tsP <- readPeakFile(tsPFile)
tsM <- readPeakFile(tsMFile)


# promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# tagMatrix <- getTagMatrix(peak, windows=promoter)
#
# data("tagMatrixList")
# tagMatrix <- tagMatrixList[[4]]

conAnn <- annotatePeak(con, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")


indAnn <- annotatePeak(ind, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db")

nonAnn <- annotatePeak(non, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

othAnn <- annotatePeak(oth, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db")


tsPAnn <- annotatePeak(tsP, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db")


tsMAnn <- annotatePeak(tsM, tssRegion=c(-3000, 3000),
                       TxDb=txdb, annoDb="org.Hs.eg.db")




# as.GRanges(conAnn)
# as.data.frame(conAnn)


write.table(x = as.data.frame(conAnn), file = conAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = as.data.frame(indAnn), file = indAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = as.data.frame(nonAnn), file = nonAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = as.data.frame(othAnn), file = othAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = as.data.frame(tsPAnn), file = tsPAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = as.data.frame(tsMAnn), file = tsMAnnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)





pdf("/kuacc/users/ualtintas20/github/Figures/CisGraph/Vers1.0/ChipSeeker.pdf", width=8.5, height=5,onefile=T)
covplot(con)
covplot(ind)
covplot(non)


plotAnnoPie(conAnn)
plotAnnoPie(indAnn)
plotAnnoPie(nonAnn)


upsetplot(conAnn)
upsetplot(indAnn)
upsetplot(nonAnn)

dev.off()
