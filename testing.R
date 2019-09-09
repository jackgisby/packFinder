source("packages.R")
source("algorithmAssessment.R")
source("devFunctions.R")

#get Arath genome
Genome <- initialise(genomeName = "Arabidopsis thaliana")

#find packs
potentialPackList <- assessPotentialPackList(subSeqs = DNAStringSet(c(1 = "CACTACAA-AAATAT",
                                                                      2 = "CACTACAA-AAATAT",
                                                                      3 = "CACTACAA-AAATAT",
                                                                      1 = "CACTACAA-AAATA",
                                                                      1 = "CACTACAA-AAA",
                                                                      1 = "CACTACAA",
                                                                      0 = "CACTACAA")),
                                              Genome = Genome,
                                              element.length = c(300, 3500),
                                              TSD.length = 3)
#filter using repeatmap


#filter using blast
