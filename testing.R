source("packages.R")

#get genome
Genome <- getGenomeDnaStringSet(genomeName = "Brassica rapa")[1:10]

#Arath
#Arabidopsis lyrata
#Arabidopsis halleri - no refseq?
#Brassica rapa

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
