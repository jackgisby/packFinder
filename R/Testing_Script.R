initialise <- function() {
  library(Biostrings)
  library(biomartr)
  library(GenomicRanges)
  library(dplyr)

  subSeq <- DNAString("CACTACAA-AAA")
  Genome <- read_genome(getGenome(db = "refseq", "Arabidopsis thaliana", path = "/Input"))
  Genome <- Genome[1:5]
}

algorithmAssessment <- function(transposonList, Genome) {
  knownCACTA <- read.csv("Input/knownCACTA.csv", sep = ";")
  knownCACTA <- knownCACTA[,1:8]
  x <- data.frame(name = Genome@ranges@NAMES)
  
  for(i in 1:length(knownCACTA$Chr)) {
    knownCACTA$chrName[i] <- as.character(x$name[knownCACTA$Chr[i]])
  }
  
  knownCACTA$Identified <- knownCACTA$start %in% transposonList$start
  
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", sum(knownCACTA$Identified)))
  print(paste0("Algorithm error rate: ", (1-(sum(knownCACTA$Identified)/length(transposonList[,1])))))
}

#initialise()
source("R/packSearch.R")

start = Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = 2, element.length = c(300, 5000))
algorithmAssessment(potentialPacks, Genome)
end = Sys.time()
print(end-start)