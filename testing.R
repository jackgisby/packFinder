source("packages.R")

#get Arath genome
Genome <- initialise(genomeName = "Arabidopsis thaliana")
subSeq <- DNAString("CACTACAA")
potentialPacks <- packSearch(subSeq, Genome, element.length = c(300, 3500), TSD.length = 3)

#find packs
potentialPackList <- assessPotentialPackList(subSeqs = DNAStringSet(c("1" = "CACTACAA-AAATAT",
                                                                      "2" = "CACTACAA-AAATAT",
                                                                      "3" = "CACTACAA-AAATAT",
                                                                      "1" = "CACTACAA-AAATA",
                                                                      "1" = "CACTACAA-AAA",
                                                                      "1" = "CACTACAA",
                                                                      "0" = "CACTACAA")),
                                              Genome = Genome,
                                              element.length = c(300, 3500),
                                              TSD.length = 3)
#filter using repeatmap
repeatMaps <- getRepeatMaps(Genome)
potentialPacks <- filterElements(potentialPacks, repeatMaps)

knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
knownCACTA <- saveReport(filter(potentialPacks, isTransposon == FALSE), subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
print(end-start)

#filter using blast
