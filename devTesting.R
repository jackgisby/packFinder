library(packFinder)

data("arabidopsisThalianaRefseq")
subSeq <- Biostrings::DNAString("CACTACAA")

for(file in 1:length(list.files("R/"))) {
  source(paste0("R/", list.files("R/")[file]))
}

packMatches <- packSearch(subSeq,
                          arabidopsisThalianaRefseq,
                          mismatch = 0,
                          elementLength = c(300, 3500),
                          tsdLength = 3)

packMatches <- packClust(packMatches,
                         arabidopsisThalianaRefseq,
                         saveFolder = "devData/",
                         vSearchPath = "D:/vsearch-2.14.1-win-x86_64/vsearch.exe")

tirClust(packMatches = packMatches,
        Genome = arabidopsisThalianaRefseq)

x <- getPackSeqs(packMatches = packMatches,
            arabidopsisThalianaRefseq)
