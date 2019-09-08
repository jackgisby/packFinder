source("R/packSearch.R")
source("R/devFunctions.R")

#pack find
Genome <- initialise()
subSeq <- DNAString("CACTACAA") #CACTACAA-AAATAT / DNAString(consensusString(knownTIRs))
max.mismatch = 0

start <- Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = max.mismatch, element.length = c(300, 3500), TSD.length = 3)
end <- Sys.time()

knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
print(end-start)

db <- blast(db="C:/Users/jackg/Documents/R/nt_db/nt/nt", type = "blastn")

#blast
start <- Sys.time()
getBlastMatches <- packBlast(potentialPacks, "nt")
end <- Sys.time()
print(end-start)
