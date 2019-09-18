source("packages.R")

genomeName <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Arabidopsis halleri", "Brassica rapa")
db <- c("refseq", "refseq", "genbank", "refseq")

#get genome
i <- 1
Genome <- getGenomeDnaStringSet(genomeName = genomeName[i], db = db[i])
sum(Genome@ranges@width)

#find packs
# assessPotentialPackList(subSeqs = DNAStringSet(c("0" = "CACTACAA")),
#                         Genome = Genome,
#                         element.length = c(300, 3500),
#                         TSD.length = 3)

Genome <- getGenomeDnaStringSet()
knownCACTA <- getArathCACTA(Genome)
CACTACAA_Data <- read.csv("Data/Output/algorithmAssessment/Full_CACTACAA_Data.csv")

#clustering
consensusSeqs <- CACTACAA_Data %>% 
  filter(Genome == unique(Genome)[[1]]) %>%
  getClusterConsensus(getOrganismHClust(., .$Genome), h = 0.69) %>% 
  print()
