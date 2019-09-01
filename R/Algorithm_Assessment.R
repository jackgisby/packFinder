#Genome <- initialise()
subSeq <- DNAString("CACTACAA-AAATAT") #CACTACAA-AAATAT
source("R/packSearch.R")

start = Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = 2, element.length = c(300, 3000), TSD.length = 3)
end = Sys.time()
identifiedCACTA <- algorithmAssessment(potentialPacks, Genome)
print(end-start)


initialise <- function() {
  "
  loads the ArAth genome and required packages for testing
  "
  
  library(Biostrings)
  library(biomartr)
  library(GenomicRanges)
  library(dplyr)

  
  Genome <- read_genome(getGenome(db = "refseq", "Arabidopsis thaliana", path = "/Input"))
  return(Genome[1:5])
}

algorithmAssessment <- function(transposonList, Genome) {
  chrNames <- data.frame(name = Genome@ranges@NAMES)
  
  knownCACTA <- read.csv("Input/knownCACTA.csv", sep = ";")[,1:8] %>% mutate(
    identified = (start %in% transposonList$start & end %in% transposonList$end)
  )
  
  identifiedCACTA <- filter(transposonList, 
                            transposonList$start %in% knownCACTA$start &
                            transposonList$end %in% knownCACTA$end) #does not consider chromosome
  
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", length(identifiedCACTA[,1])))
  print(paste0("Algorithm error rate: ", (1-(length(identifiedCACTA[,1])/length(transposonList[,1])))))
  
  return(identifiedCACTA)
}

