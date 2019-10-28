# this script was used to generate the file found at:
# inst/extdata/packMatches.blast6out

library(packFinder)

data("arabidopsisThalianaRefseq")
data("packMatches")

packClusts <- packClust(
    packMatches, 
    arabidopsisThalianaRefseq,
    saveFolder = system.file("extdata", package = "packFinder")
)
