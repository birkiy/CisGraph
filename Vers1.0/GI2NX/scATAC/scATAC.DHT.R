
.libPaths( c( "/scratch/users/ualtintas20/apps/R/lib64", "/scratch/kuacc/apps/R/3.6.1/lib64/R/library", .libPaths() ) )


suppressMessages(library(GenomicInteractions))

library(rtracklayer)
library(InteractionSet)


# home = "/home/birkiy/github/CisGraph/Vers1.0"
home = "/kuacc/users/ualtintas20/github/Data/CisGraph/Vers1.0/InteractionBedPe"

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




InteractionFile = paste(home, "LNCaP_DHT_2000_CiceroConns.05Filter.bedpe", sep="/")
Conn.rep1 = makeGInteractionsFromGRangesPairs(import(InteractionFile))
# = makeGenomicInteractionsFromFile(InteractionFile,
# type="chiapet.tool", experiment_name="LNCaP.DHT", description="LNCaP.DHT.Cicero")

conFile = paste(home, "Data/Regions/cons-arbs.bed", sep="/")
conBed = bed_to_granges(conFile)

indFile = paste(home, "Data/Regions/ind-arbs.bed", sep="/")
indBed = bed_to_granges(indFile)

nonFile = paste(home, "Data/Regions/Non-Active-ARBS.bed", sep="/")
nonBed = bed_to_granges(nonFile)


proFile = paste(home, "Data/Regions/promoters_ann_5kb.bed", sep="/")
proBed = bed_to_granges(proFile)


print(
"***
Annotation Interaction Starts
***"
)

annotation.features = list(pro = proBed, con = conBed, ind = indBed, non = nonBed)

annotateInteractions(Conn.rep1, annotation.features)

print(
"***
Remove Distal
***"
)



Conn.rep1 <- Conn.rep1[!(anchorOne(Conn.rep1)$node.class == "distal")]

Conn <- Conn.rep1[!(anchorTwo(Conn.rep1)$node.class == "distal")]


print(
  categoriseInteractions(Conn)
)


print(
"***
Start Converting the GI object that compatible with Cis-Graph-NetworkX
***"
)

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


df_all = as.data.frame(Conn[isInteractionType(Conn, c("pro", "con", "ind", "non"),
                                              c("pro", "con", "ind", "non"))])

print(
"***
Multiple extraction
***"
)



# pro

df_all <- df_all %>% separate_rows(pro.id1,sep = "\"")
df_all <- df_all %>% separate_rows(pro.id2,sep = "\"")
df_all <- df_all[df_all$pro.id1!="c(",]
df_all <- df_all[df_all$pro.id2!="c(",]
df_all <- df_all[df_all$pro.id1!=", ",]
df_all <- df_all[df_all$pro.id2!=", ",]
df_all <- df_all[df_all$pro.id1!=")",]
df_all <- df_all[df_all$pro.id2!=")",]
df_all <- df_all[df_all$pro.id1!="",]
df_all <- df_all[df_all$pro.id2!="",]

df_all$pro.id1 <- paste("pro",df_all$pro.id1,sep=".")
df_all$pro.id2 <- paste("pro",df_all$pro.id2,sep=".")



# con

df_all <- df_all %>% separate_rows(con.id1,sep = "\"")
df_all <- df_all %>% separate_rows(con.id2,sep = "\"")
df_all <- df_all[df_all$con.id1!="c(",]
df_all <- df_all[df_all$con.id2!="c(",]
df_all <- df_all[df_all$con.id1!="",]
df_all <- df_all[df_all$con.id2!="",]
df_all <- df_all[df_all$con.id1!=", ",]
df_all <- df_all[df_all$con.id2!=", ",]
df_all <- df_all[df_all$con.id1!=")",]
df_all <- df_all[df_all$con.id2!=")",]

df_all$con.id1 <- paste("con",df_all$con.id1,sep=".")
df_all$con.id2 <- paste("con",df_all$con.id2,sep=".")


# ind

df_all <- df_all %>% separate_rows(ind.id1,sep = "\"")
df_all <- df_all %>% separate_rows(ind.id2,sep = "\"")
df_all <- df_all[df_all$ind.id1!="c(",]
df_all <- df_all[df_all$ind.id2!="c(",]
df_all <- df_all[df_all$ind.id1!="",]
df_all <- df_all[df_all$ind.id2!="",]
df_all <- df_all[df_all$ind.id1!=", ",]
df_all <- df_all[df_all$ind.id2!=", ",]
df_all <- df_all[df_all$ind.id1!=")",]
df_all <- df_all[df_all$ind.id2!=")",]

df_all$ind.id1 <- paste("ind",df_all$ind.id1,sep=".")
df_all$ind.id2 <- paste("ind",df_all$ind.id2,sep=".")


# non

df_all <- df_all %>% separate_rows(non.id1,sep = "\"")
df_all <- df_all %>% separate_rows(non.id2,sep = "\"")
df_all <- df_all[df_all$non.id1!="c(",]
df_all <- df_all[df_all$non.id2!="c(",]
df_all <- df_all[df_all$non.id1!="",]
df_all <- df_all[df_all$non.id2!="",]
df_all <- df_all[df_all$non.id1!=", ",]
df_all <- df_all[df_all$non.id2!=", ",]
df_all <- df_all[df_all$non.id1!=")",]
df_all <- df_all[df_all$non.id2!=")",]

df_all$non.id1 <- paste("non",df_all$non.id1,sep=".")
df_all$non.id2 <- paste("non",df_all$non.id2,sep=".")



print(
"***
Concatanate different classes to seperate later
***"
)

a = paste(df_all$pro.id1, df_all$con.id1, df_all$ind.id1, df_all$non.id1, sep=",")


b = paste(df_all$pro.id2, df_all$con.id2, df_all$ind.id2, df_all$non.id2, sep=",")



print(
"***
Create df3 and seperate rows to remove NAs later
***"
)

df3 <- data.frame(aNodes = a, bNodes = b)

df3 <- df3 %>% separate_rows(aNodes,sep = ",")
df3 <- df3 %>% separate_rows(bNodes,sep = ",")


print(
"***
Remove rows if they have NA
***"
)

df4 = df3[(df3$aNodes %in% "con.NA" + df3$bNodes %in% "con.NA"
	+ df3$aNodes %in% "con.NA" + df3$bNodes %in% "con.NA"
  + df3$aNodes %in% "ind.NA" + df3$bNodes %in% "ind.NA"
  + df3$aNodes %in% "non.NA" + df3$bNodes %in% "non.NA"
  + df3$aNodes %in% "pro.NA" + df3$bNodes %in% "pro.NA") == 0,]

df4 = df4[!(df4[,2] == ""),]



print(
"***
Concatanate aNodes and bNodes and aggregate
***"
)

df4 = data.frame(abNodes = paste(df4$aNodes, df4$bNodes, sep=","))


#
# df5 <- cbind(
#   aggregate(df4$weight, by=list(abNodes=df4$abNodes), sum)[,c(1,2)],
#   aggregate(df4$fdr, by=list(abNodes=df4$abNodes), mean)[,2]
# )
#

print(
"***
Write table
***"
)



write.table(df4,file = "Data/GIs/GI.G.scATAC.DHT.txt",
            sep="\t", quote = F, col.names = F, row.names = F)
