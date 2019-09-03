source("R/packSearch.R")
source("R/devFunctions.R")

#Genome <- initialise()
subSeq <- DNAString("CACTACAA") #CACTACAA-AAATAT / DNAString(consensusString(knownTIRs))
max.mismatch = 0

start <- Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = max.mismatch, element.length = c(300, 3500), TSD.length = 3)
end <- Sys.time()

knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
print(end-start)
