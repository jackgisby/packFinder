source("R/devFunctions.R")
subSeq <- DNAString("CACTACAA-AAATAT") #CACTACAA-AAATAT

knownCACTA <- getArAthCACTA(Genome)
knownTIRs <- DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR))

badMatches <- which(!assessSubSeq(subSeq, knownTIRs, mismatch = 7))
badMatches <- knownTIRs[badMatches]
print(length(badMatches))
print(badMatches)