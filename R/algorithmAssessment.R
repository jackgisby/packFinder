source("R/packSearch.R")
source("R/devFunctions.R")

#Genome <- initialise()
subSeq <- DNAString("CACTACAA-AAATAT") #CACTACAA-AAATAT / DNAString(consensusString(knownTIRs))


start <- Sys.time()
potentialPacks <- packSearch(subSeq, Genome, mismatch = 2, element.length = c(300, 3500), TSD.length = 3)
end <- Sys.time()

identifiedCACTA <- algorithmAssessment(potentialPacks, Genome)
knownCACTA <- getArAthCACTA(Genome)
knownTIRs <- DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR))
badCACTA <- getBadMatches(knownCACTA, subSeq, 2)
print(end-start)
