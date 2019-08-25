if(initialise == TRUE) {
  library(Biostrings)
  library(biomartr)
  library(GenomicRanges)
  library(dplyr)

  subSeq <- DNAString("CACTACAA-AAA")
  Genome <- read_genome(getGenome(db = "refseq", "Arabidopsis thaliana", path = "/Input"))
  Genome <- Genome[1:5]
}

initialise <- FALSE

source("R/packSearch.R")

start = Sys.time()
forwardMatches <- packSearch(subSeq, Genome, mismatch = 2, element.length = c(300, 5000))
#x <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, 1, c(300, 5000), 3)
end = Sys.time()
print(end-start)