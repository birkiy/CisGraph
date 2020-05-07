

suppressMessages(library(GenomicInteractions))


print(
"***
Bed reader
***"
)

bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)

  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }


  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  }
    return(gr)

}



print(
"***
G, the basal graph
***"
)

home = "/home/birkiy/github/CisGraph"

InteractionFile = paste(home, "Data/InteractionBedPe/VCaP_AR_ChIA-PET.bedpe", sep="/")

VCaP.rep1 = makeGenomicInteractionsFromFile(InteractionFile,
	type="chiapet.tool", experiment_name="VCaP", description="ChiA-PET VCaP")

enhFile = paste(home, "Data/Regions/EnhancerDomain.bed", sep="/")
enhBed = bed_to_granges(enhFile)

enPFile = paste(home, "Data/Regions/EPDomain.bed", sep="/")
enPBed = bed_to_granges(enPFile)

genFile = paste(home, "Data/Regions/PromoterDomain.bed", sep="/")
genBed = bed_to_granges(genFile)


print(
"***
Annotation Interaction Starts
***"
)

annotation.features = list(gen = genBed, enP = enPBed, enh = enhBed)

annotateInteractions(VCaP.rep1, annotation.features)

print(
"***
Remove Distal
***"
)



VCaP.rep1 <- VCaP.rep1[!(anchorOne(VCaP.rep1)$node.class == "distal")]

VCaP <- VCaP.rep1[!(anchorTwo(VCaP.rep1)$node.class == "distal")]


print(
  categoriseInteractions(VCaP)
)


print(
"***
Start Converting the GI object that compatible with Cis-Graph-NetworkX
***"
)

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


df_all = as.data.frame(VCaP[isInteractionType(VCaP, c("gen", "enP", "enh"),
                                              c("gen", "enP", "enh"))])

print(
"***
Multiple extraction
***"
)



# gen

df_all <- df_all %>% separate_rows(gen.id1,sep = "\"")
df_all <- df_all %>% separate_rows(gen.id2,sep = "\"")
df_all <- df_all[df_all$gen.id1!="c(",]
df_all <- df_all[df_all$gen.id2!="c(",]
df_all <- df_all[df_all$gen.id1!=", ",]
df_all <- df_all[df_all$gen.id2!=", ",]
df_all <- df_all[df_all$gen.id1!=")",]
df_all <- df_all[df_all$gen.id2!=")",]
df_all <- df_all[df_all$gen.id1!="",]
df_all <- df_all[df_all$gen.id2!="",]

df_all$gen.id1 <- paste("gen",df_all$gen.id1,sep=".")
df_all$gen.id2 <- paste("gen",df_all$gen.id2,sep=".")



# enP

df_all <- df_all %>% separate_rows(enP.id1,sep = "\"")
df_all <- df_all %>% separate_rows(enP.id2,sep = "\"")
df_all <- df_all[df_all$enP.id1!="c(",]
df_all <- df_all[df_all$enP.id2!="c(",]
df_all <- df_all[df_all$enP.id1!="",]
df_all <- df_all[df_all$enP.id2!="",]
df_all <- df_all[df_all$enP.id1!=", ",]
df_all <- df_all[df_all$enP.id2!=", ",]
df_all <- df_all[df_all$enP.id1!=")",]
df_all <- df_all[df_all$enP.id2!=")",]

df_all$enP.id1 <- paste("enP",df_all$enP.id1,sep=".")
df_all$enP.id2 <- paste("enP",df_all$enP.id2,sep=".")


# enh

df_all <- df_all %>% separate_rows(enh.id1,sep = "\"")
df_all <- df_all %>% separate_rows(enh.id2,sep = "\"")
df_all <- df_all[df_all$enh.id1!="c(",]
df_all <- df_all[df_all$enh.id2!="c(",]
df_all <- df_all[df_all$enh.id1!="",]
df_all <- df_all[df_all$enh.id2!="",]
df_all <- df_all[df_all$enh.id1!=", ",]
df_all <- df_all[df_all$enh.id2!=", ",]
df_all <- df_all[df_all$enh.id1!=")",]
df_all <- df_all[df_all$enh.id2!=")",]

df_all$enh.id1 <- paste("enh",df_all$enh.id1,sep=".")
df_all$enh.id2 <- paste("enh",df_all$enh.id2,sep=".")



print(
"***
Concatanate different classes to seperate later
***"
)


a = paste(df_all$gen.id1, df_all$enP.id1, df_all$enh.id1,sep=",")


b = paste(df_all$gen.id2, df_all$enP.id2, df_all$enh.id2,sep=",")



print(
"***
Create df3 and seperate rows to remove NAs later
***"
)

df3 <- data.frame(aNodes = a, bNodes = b, weight = df_all$counts, fdr = df_all$fdr)

df3 <- df3 %>% separate_rows(aNodes,sep = ",")
df3 <- df3 %>% separate_rows(bNodes,sep = ",")


print(
"***
Remove rows if they have NA
***"
)

df4 = df3[(df3$aNodes %in% "gen.NA" + df3$bNodes %in% "gen.NA"
	+ df3$aNodes %in% "enP.NA" + df3$bNodes %in% "enP.NA"
  + df3$aNodes %in% "enh.NA" + df3$bNodes %in% "enh.NA") == 0,]

df4 = df4[!(df4[,2] == ""),]



print(
"***
Concatanate aNodes and bNodes and aggregate
***"
)

df4 = data.frame(abNodes = paste(df4$aNodes, df4$bNodes, sep=","), weight = df4$weight, fdr = df4$fdr)



df5 <- cbind(
  aggregate(df4$weight, by=list(abNodes=df4$abNodes), sum)[,c(1,2)],
  aggregate(df4$fdr, by=list(abNodes=df4$abNodes), mean)[,2]
)


print(
"***
Write table
***"
)



write.table(df5,file = "Data/GIs/GI.G.txt",
            sep="\t", quote = F, col.names = F, row.names = F)
