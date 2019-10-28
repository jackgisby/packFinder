# this script was used to generate the file found at:
# inst/extdata/packMatches.fasta

library(packFinder)

data("arabidopsisThalianaRefseq")
data("packMatches")

packsToFasta(packMatches, system.file(
    "extdata", 
    "packMatches.fasta", 
    package = "packFinder"
), arabidopsisThalianaRefseq)
