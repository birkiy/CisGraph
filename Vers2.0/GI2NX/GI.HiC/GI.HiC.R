


GIobj <- function(home, interactionFile, bedFile, nodeClass, outFile){
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
  The basal graph
  ***"
  )

  InteractionFile = paste(home, interactionFile, sep="/")

  VCaP.rep1 = makeGenomicInteractionsFromFile(InteractionFile,
  	type="chiapet.tool", experiment_name="VCaP", description="ChiA-PET VCaP")


  bedFile = paste(home, bedFile, sep="/")
  nodeBed = bed_to_granges(bedFile)


  print(
  "***
  Annotation Interaction Starts
  ***"
  )

  annotation.features = list(node = nodeBed
                             )

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



  df_all = as.data.frame(VCaP[isInteractionType(VCaP, c("node"),
                                                c("node"))])



  print(
  "***
  Multiple extraction
  ***"
  )



  # tad

  df_all <- df_all %>% separate_rows(node.id1,sep = "\"")
  df_all <- df_all %>% separate_rows(node.id2,sep = "\"")
  df_all <- df_all[df_all$node.id1!="c(",]
  df_all <- df_all[df_all$node.id2!="c(",]
  df_all <- df_all[df_all$node.id1!=", ",]
  df_all <- df_all[df_all$node.id2!=", ",]
  df_all <- df_all[df_all$node.id1!=")",]
  df_all <- df_all[df_all$node.id2!=")",]
  df_all <- df_all[df_all$node.id1!="",]
  df_all <- df_all[df_all$node.id2!="",]

  df_all$node.id1 <- paste(nodeClass,df_all$node.id1,sep=".")
  df_all$node.id2 <- paste(nodeClass,df_all$node.id2,sep=".")





  print(
  "***
  Concatanate aNodes and bNodes and aggregate
  ***"
  )


  a = df_all$node.id1


  b = df_all$node.id2



  df3 <- data.frame(aNodes = a, bNodes = b, weight = df_all$counts, fdr = df_all$fdr)
  df4 = df3[(df3$aNodes %in% paste(nodeClass, "NA", sep=".") + df3$bNodes %in% paste(nodeClass, "NA", sep=".")) == 0,]
  df4 = df4[!(df4[,2] == ""),]



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

  write.table(df5,file = paste(home, outFile, sep="/"), sep="\t", quote = F, col.names = F, row.names = F)


}



.libPaths(c("/home/ualtintas/anaconda3/envs/CisGraph/lib/R/library", .libPaths()))


dataPwd = "/home/ualtintas/github/Data/Vers2.0"
interactionFile = "InteractionBedPe/FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bedpe"

print(
  "****
  G starts
  ****")
GIobj(home=dataPwd, interactionFile=interactionFile, bedFile="NodeBeds/creCtr.bed", outFile="GIs/GI.G.txt", nodeClass="cre")



# print(
#   "****
#   T starts
#   ****")
# GIobj(home=home, interactionFile=interactionFile, bedFile="Data/Regions/TAD.C42B.bed", outFile="Data/GIs/GI.T.txt", nodeClass="tad")
#
#
#
# print(
#   "****
#   C starts
#   ****")
# GIobj(home=home, interactionFile=interactionFile, bedFile="Data/Regions/hg19.Chr.bed", outFile="Data/GIs/GI.C.txt", nodeClass="chr")
