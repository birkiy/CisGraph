

# library("tximport")
# library("readr")
# library("tximportData")
# dir <- system.file("extdata", package="tximportData")
# samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# samples$condition <- factor(rep(c("A","B"),each=3))
# rownames(samples) <- samples$run
# samples[,c("pop","center","run","condition")]



totalMatrix = read.csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/countTableRaw.txt", sep="\t")

samples = c(
  "GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5",
  "GR.12h.dex.rep1", "GR.12h.dex.rep2", "GR.12h.dex.rep3", "GR.12h.dex.rep4", "GR.12h.dex.rep5",
  "GR.1h.dex.rep1", "GR.1h.dex.rep2", "GR.1h.dex.rep3", "GR.1h.dex.rep4", "GR.1h.dex.rep5",
  "GR.4h.dex.rep1", "GR.4h.dex.rep2", "GR.4h.dex.rep3", "GR.4h.dex.rep4", "GR.4h.dex.rep5",
  "GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5",
  "input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9")



mat = as.matrix(totalMatrix[,4:40])

rownames(mat) = apply(totalMatrix[,1:3], 1, paste, collapse="-")
colnames(mat) = samples

totalColData = data.frame(condition=c(rep("0h", 5), rep("12h", 5), rep("1h", 5), rep("4h", 5), rep("8h", 5), rep("input", 12)), lib=c(rep("starr", 25), rep("input", 12)))
rownames(totalColData) = samples
totalColData$condition <- factor(totalColData$condition)
totalColData$lib <- factor(totalColData$lib)


library("DESeq2") # => Note to activate r-environment

dds <- DESeqDataSetFromMatrix(countData = mat[,1:25],
                              colData = totalColData[1:25,],
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "0h")
dds <- DESeq(dds)
res <- results(dds)



res1hDex = results(dds, contrast=c("condition", "1h","0h"), alpha=0.05)
res4hDex = results(dds, contrast=c("condition", "4h","0h"), alpha=0.05)
res8hDex = results(dds, contrast=c("condition", "8h","0h"), alpha=0.05)
res12hDex = results(dds, contrast=c("condition", "12h","0h"), alpha=0.05)




pdf("MAplotsGRdex.pdf", width = 6.6, height = 6)

par(mfrow=c(2,2), mar=c(4,4,2,1))
plotMA(res1hDex, ylim=c(-2,2), main="1h (padj<0.05)")
plotMA(res4hDex, ylim=c(-2,2), main="4h (padj<0.05)")
plotMA(res8hDex, ylim=c(-2,2), main="8h (padj<0.05)")
plotMA(res12hDex, ylim=c(-2,2), main="12h (padj<0.05)")

dev.off()



write.table(as.data.frame(res1hDex), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.1h.dex.tsv")

write.table(as.data.frame(res4hDex), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.4h.dex.tsv")

write.table(as.data.frame(res8hDex), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.8h.dex.tsv")

write.table(as.data.frame(res12hDex), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.12h.dex.tsv")



################################333333
dds <- DESeqDataSetFromMatrix(countData = mat[,1:25],
                              colData = totalColData[1:25,],
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "0h")
dds <- DESeq(dds)
res <- results(dds)




res8hDex = results(dds, contrast=c("lib", "starr","input"), alpha=0.05)

# 0h
idx0 = c(1:5, 26:37)
mat0 = mat[,idx0]
colData0 = data.frame(condition=c(rep("0h", 5), rep("input", 12)))
rownames(colData0) = samples[idx0]
colData0$condition <- factor(colData0$condition)


dds <- DESeqDataSetFromMatrix(countData = mat0,
                              colData = colData0,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("0h","input"))
dds$condition <- relevel(dds$condition, ref = "input")
# collapseReplicates => can be used to merge libraries
dds <- DESeq(dds)
res <- results(dds)
res005 <- results(dds, alpha=0.05)


write.table(as.data.frame(res005), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.0h.lib.tsv")

#######################################3

# 8h
idx8 = c(21:25, 26:37)
mat8 = mat[,idx8]
colData8 = data.frame(condition=c(rep("8h", 5), rep("input", 12)))
rownames(colData8) = samples[idx8]
colData8$condition <- factor(colData8$condition)


dds <- DESeqDataSetFromMatrix(countData = mat8,
                              colData = colData8,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("8h","input"))
dds$condition <- relevel(dds$condition, ref = "input")
# collapseReplicates => can be used to merge libraries
dds <- DESeq(dds)
res <- results(dds)
res805 <- results(dds, alpha=0.05)


write.table(as.data.frame(res805), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.8h.lib.tsv")




#####################################3




# 0h
idx0 = c(1:5, 21:25, 26:37)
mat0 = mat[,idx0]
colData0 = data.frame(condition=c(rep("starr", 10), rep("input", 12)))
rownames(colData0) = samples[idx0]
colData0$condition <- factor(colData0$condition)


dds <- DESeqDataSetFromMatrix(countData = mat0,
                              colData = colData0,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("starr","input"))
dds$condition <- relevel(dds$condition, ref = "input")
# collapseReplicates => can be used to merge libraries
dds <- DESeq(dds)
res <- results(dds)
res005 <- results(dds, alpha=0.05)


write.table(as.data.frame(res005), sep="\t",  quote = FALSE,
          file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.starr.lib.tsv")






# # 0h
# idx0 = c(1:5, 26:37)
# mat0 = mat[,idx0]
# colData0 = data.frame(condition=c(rep("0h", 5), rep("input", 12)))
# rownames(colData0) = samples[idx0]
# colData0$condition <- factor(colData0$condition)
#
# # 1h
# idx1 = c(11:15, 26:37)
# mat1 = mat[,idx1]
# colData1 = data.frame(condition=c(rep("1h", 5), rep("input", 12)))
# rownames(colData1) = samples[idx1]
# colData1$condition <- factor(colData1$condition)
#
# # 4h
# idx4 = c(16:20, 26:37)
# mat4 = mat[,idx4]
# colData4 = data.frame(condition=c(rep("4h", 5), rep("input", 12)))
# rownames(colData4) = samples[idx4]
# colData4$condition <- factor(colData4$condition)
#
# # 8h
# idx8 = c(21:25, 26:37)
# mat8 = mat[,idx8]
# colData8 = data.frame(condition=c(rep("8h", 5), rep("input", 12)))
# rownames(colData8) = samples[idx8]
# colData8$condition <- factor(colData8$condition)
#
# # 12h
# idx12 = c(6:10, 26:37)
# mat12 = mat[,idx12]
# colData12 = data.frame(condition=c(rep("12h", 5), rep("input", 12)))
# rownames(colData12) = samples[idx12]
# colData12$condition <- factor(colData12$condition)
#
#
#
#
#
# dds <- DESeqDataSetFromMatrix(countData = mat0,
#                               colData = colData0,
#                               design = ~ condition)
# dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("0h","input"))
# dds$condition <- relevel(dds$condition, ref = "input")
# # collapseReplicates => can be used to merge libraries
# dds <- DESeq(dds)
# res <- results(dds)
# res005 <- results(dds, alpha=0.05)
#
#
# dds <- DESeqDataSetFromMatrix(countData = mat1,
#                               colData = colData1,
#                               design = ~ condition)
# dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("1h","input"))
# dds$condition <- relevel(dds$condition, ref = "input")
# # collapseReplicates => can be used to merge libraries
# dds <- DESeq(dds)
# res <- results(dds)
# res105 <- results(dds, alpha=0.05)
#
#
# dds <- DESeqDataSetFromMatrix(countData = mat4,
#                               colData = colData4,
#                               design = ~ condition)
# dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("4h","input"))
# dds$condition <- relevel(dds$condition, ref = "input")
# # collapseReplicates => can be used to merge libraries
# dds <- DESeq(dds)
# res <- results(dds)
# res405 <- results(dds, alpha=0.05)
#
#
# dds <- DESeqDataSetFromMatrix(countData = mat8,
#                               colData = colData8,
#                               design = ~ condition)
# dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("8h","input"))
# dds$condition <- relevel(dds$condition, ref = "input")
# # collapseReplicates => can be used to merge libraries
# dds <- DESeq(dds)
# res <- results(dds)
# res805 <- results(dds, alpha=0.05)
#
#
# dds <- DESeqDataSetFromMatrix(countData = mat12,
#                               colData = colData12,
#                               design = ~ condition)
# dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("12h","input"))
# dds$condition <- relevel(dds$condition, ref = "input")
# # collapseReplicates => can be used to merge libraries
# dds <- DESeq(dds)
# res <- results(dds)
# res1205 <- results(dds, alpha=0.05)
#
#
# pdf("MAplotsGRstarr.pdf", width = 10, height = 6)
#
# par(mfrow=c(2,3), mar=c(4,4,2,1))
# plotMA(res005, ylim=c(-2,2), main="0h (padj<0.05)")
# plotMA(res105, ylim=c(-2,2), main="1h (padj<0.05)")
# plotMA(res405, ylim=c(-2,2), main="4h (padj<0.05)")
# plotMA(res805, ylim=c(-2,2), main="8h (padj<0.05)")
# plotMA(res1205, ylim=c(-2,2), main="12h (padj<0.05)")
#
# dev.off()
#
# write.table(as.data.frame(res005), sep="\t",  quote = FALSE,
#           file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.0h.deg.tsv")
#
#
#
# write.table(as.data.frame(res105), sep="\t",  quote = FALSE,
#           file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.1h.deg.tsv")
#
#
# write.table(as.data.frame(res405), sep="\t",  quote = FALSE,
#           file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.4h.deg.tsv")
#
#
# write.table(as.data.frame(res805), sep="\t",  quote = FALSE,
#           file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.8h.deg.tsv")
#
#
# write.table(as.data.frame(res1205), sep="\t",  quote = FALSE,
#           file="/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.12h.deg.tsv")
