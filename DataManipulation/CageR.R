



library(CAGEr)

inputFiles = c("/home/ualtintas/LncapExo/CAGE/cage.dht.bam",
               "/home/ualtintas/LncapExo/CAGE/cage.etoh.bam")


inputFiles = list.files( system.file("extdata", package = "CAGEr")
                       , "ctss$"
                       , full.names = TRUE)


ce <- CAGEexp(genomeName     = "BSgenome.Drerio.UCSC.danRer7" ,
              inputFiles     = inputFiles,
              inputFilesType = "bam" ,
              sampleLabels   = sub( ".chr17.ctss", "", basename(inputFiles))
)
