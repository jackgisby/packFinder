library(packFinder)

data("arabidopsisThalianaRefseq")
subSeq <- Biostrings::DNAString("CACTACAA")

packMatches <- packSearch(subSeq,
                          arabidopsisThalianaRefseq,
                          mismatch = 0,
                          elementLength = c(300, 3500),
                          tsdLength = 3)

packMatches <- packClust(packMatches,
                         saveFolder = "devData/",
                         vSearchPath = "D:/vsearch-2.14.1-win-x86_64/vsearch.exe")
