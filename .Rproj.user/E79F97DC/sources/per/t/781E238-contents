source("packages.R")

genomeName <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Arabidopsis halleri", "Brassica rapa")
db <- c("refseq", "refseq", "genbank", "refseq")

#get genome
i <- 4
Genome <- getGenomeDnaStringSet(genomeName = genomeName[i], db = db[i])


#find packs
assessPotentialPackList(subSeqs = DNAStringSet(c("1" = "CACTACAA-AAATAT",
                                                 "2" = "CACTACAA-AAATAT",
                                                 "1" = "CACTACAA-AAATA",
                                                 "1" = "CACTACAA-AAA",
                                                 "0" = "CACTACAA")),
                        Genome = Genome,
                        element.length = c(300, 3500),
                        TSD.length = 3)

#clustering
potentialPacks <- read.csv("Data/Output/algorithmAssessment/potentialPacks.csv") %>%
  filter(stringID == 5) %>%
  mutate(forward_TIR = as.character(forward_TIR)) %>%
  mutate(reverse_TIR = as.character(reverse_TIR))

getTirClusters(potentialPacks, model = "TN93")
