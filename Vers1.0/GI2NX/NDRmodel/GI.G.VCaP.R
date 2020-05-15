

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

dhsFile = paste(home, "Data/Regions/creDHS.bed", sep="/")
dhsBed = bed_to_granges(dhsFile)


print(
"***
Annotation Interaction Starts
***"
)

annotation.features = list(cre = dhsBed)

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


df_all = as.data.frame(VCaP[isInteractionType(VCaP, c("cre"),
                                              c("cre"))])

print(
"***
Multiple extraction
***"
)




# tad

df_all <- df_all %>% separate_rows(cre.id1,sep = "\"")
df_all <- df_all %>% separate_rows(cre.id2,sep = "\"")
df_all <- df_all[df_all$cre.id1!="c(",]
df_all <- df_all[df_all$cre.id2!="c(",]
df_all <- df_all[df_all$cre.id1!=", ",]
df_all <- df_all[df_all$cre.id2!=", ",]
df_all <- df_all[df_all$cre.id1!=")",]
df_all <- df_all[df_all$cre.id2!=")",]
df_all <- df_all[df_all$cre.id1!="",]
df_all <- df_all[df_all$cre.id2!="",]

df_all$cre.id1 <- paste("cre",df_all$cre.id1,sep=".")
df_all$cre.id2 <- paste("cre",df_all$cre.id2,sep=".")





print(
"***
Concatanate aNodes and bNodes and aggregate
***"
)


df3 <- data.frame(abNodes = paste(df_all$cre.id1,df_all$cre.id2, sep=","), weight = df_all$counts, fdr = df_all$fdr)



df3 <- cbind(
  aggregate(df3$weight, by=list(abNodes=df3$abNodes), sum)[,c(1,2)],
  aggregate(df3$fdr, by=list(abNodes=df3$abNodes), mean)[,2]
)



print(
"***
Write table
***"
)

write.table(df3,file = "Data/GIs/dhs/GI.G.txt",
            sep="\t", quote = F, col.names = F, row.names = F)
