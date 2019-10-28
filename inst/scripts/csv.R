# this script was used to generate the file found at:
# inst/extdata/packMatches.csv

library(packFinder)

data("packMatches")

packsToCsv(packMatches, system.file(
    "extdata", 
    "packMatches.csv", 
    package = "packFinder"
))
