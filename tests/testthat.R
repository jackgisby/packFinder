library(testthat)
library(packFinder)

data("arabidopsisThalianaRefseq")
subSeq <- Biostrings::DNAString("CACTACAA")

packMatches <- packSearch(subSeq, Genome, mismatch = 0, elementLength = c(300, 3500), tsdLength = 3)

packClusts <- packClust(packMatches,
                        saveFolder = "devTestOutput/",
                        vSearchPath = "D:/vsearch-2.14.1-win-x86_64/vsearch.exe")

consensusSeqs <- tirClust(packClusts,
                          tirLength = 25,
                          plotSavePath = "devTestOutput/tirRelationships.png")
