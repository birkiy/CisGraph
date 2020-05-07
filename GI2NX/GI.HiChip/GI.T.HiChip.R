

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
T, the basal graph
***"
)

home = "/kuacc/users/ualtintas20/CisGraph"

InteractionFile = paste(home, "Data/InteractionBedPe/FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bedpe", sep="/")

HiChip.rep1 = makeGenomicInteractionsFromFile(InteractionFile,
	type="HiChip", experiment_name="HiChip", description="HiChip LNCaP")


tadFile = paste(home, "Data/Regions/TAD.C42B.bed", sep="/")
tadBed = bed_to_granges(tadFile)


print(
"***
Annotation Interaction Starts
***"
)

annotation.features = list(tad = tadBed
                           )

annotateInteractions(HiChip.rep1, annotation.features)


print(
"***
Remove Distal
***"
)



HiChip.rep1 <- HiChip.rep1[!(anchorOne(HiChip.rep1)$node.class == "distal")]

HiChip <- HiChip.rep1[!(anchorTwo(HiChip.rep1)$node.class == "distal")]


print(
  categoriseInteractions(HiChip)
)


print(
"***
Start Converting the GI object that compatible with Cis-Graph-NetworkX
***"
)

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))



df_all = as.data.frame(HiChip[isInteractionType(HiChip, c("tad"),
                                              c("tad"))])



print(
"***
Multiple extraction
***"
)



# tad

df_all <- df_all %>% separate_rows(tad.id1,sep = "\"")
df_all <- df_all %>% separate_rows(tad.id2,sep = "\"")
df_all <- df_all[df_all$tad.id1!="c(",]
df_all <- df_all[df_all$tad.id2!="c(",]
df_all <- df_all[df_all$tad.id1!=", ",]
df_all <- df_all[df_all$tad.id2!=", ",]
df_all <- df_all[df_all$tad.id1!=")",]
df_all <- df_all[df_all$tad.id2!=")",]
df_all <- df_all[df_all$tad.id1!="",]
df_all <- df_all[df_all$tad.id2!="",]

df_all$tad.id1 <- paste("tad",df_all$tad.id1,sep=".")
df_all$tad.id2 <- paste("tad",df_all$tad.id2,sep=".")





print(
"***
Concatanate aNodes and bNodes and aggregate
***"
)


df3 <- data.frame(abNodes = paste(df_all$tad.id1,df_all$tad.id2, sep=","), weight = df_all$counts, fdr = df_all$fdr)



df3 <- cbind(
  aggregate(df3$weight, by=list(abNodes=df3$abNodes), sum)[,c(1,2)],
  aggregate(df3$fdr, by=list(abNodes=df3$abNodes), mean)[,2]
)



print(
"***
Write table
***"
)

write.table(df3,file = "Data/GIs/HiChip/GI.T.txt",
            sep="\t", quote = F, col.names = F, row.names = F)
