source("packages.R")

Genome <- getGenomeDnaStringSet()
subSeq <- DNAString("CACTACAA")
packMatches <- packSearch(subSeq, Genome, elementLength = c(300,3500), tsdLength = 3)

packMatches <- knownCacta(Genome)
packClusts <- packClust(packMatches, Genome, identity = 0.5)
packClusts$cluster == correctClusts


source("packages.R")
consensusSeqs <- tirClust(packClusts, tirLength = 100)
consensusSeqs
