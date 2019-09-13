source("packages.R")
# 
# genomeName <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Arabidopsis halleri", "Brassica rapa")
# db <- c("refseq", "refseq", "genbank", "refseq")
# 
# #get genome
# i <- 4
# Genome <- getGenomeDnaStringSet(genomeName = genomeName[i], db = db[i])
# 
# 
# #find packs
# assessPotentialPackList(subSeqs = DNAStringSet(c("1" = "CACTACAA-AAATAT",
#                                                  "2" = "CACTACAA-AAATAT",
#                                                  "1" = "CACTACAA-AAATA",
#                                                  "1" = "CACTACAA-AAA",
#                                                  "0" = "CACTACAA")),
#                         Genome = Genome,
#                         element.length = c(300, 3500),
#                         TSD.length = 3)

Genome <- getGenomeDnaStringSet()
knownCACTA <- getArathCACTA(Genome)

#get previously generated data
getFiles <- function() {
  files <- vector("list", length(list.files("Data/Output/algorithmAssessment/")))
  
  for(file in 1:length(list.files("Data/Output/algorithmAssessment/"))) {
    files[[file]] <- read.csv(paste0("Data/Output/algorithmAssessment/",
                                     list.files("Data/Output/algorithmAssessment/")[file],
                                     "/potentialPacks.csv"))
    
    files[[file]]$Genome <- list.files("Data/Output/algorithmAssessment/")[file]
    
  }
  
  rbind(files[[1]],
        files[[2]],
        files[[3]],
        files[[4]]) %>%
    filter(TSD != "NNN") %>%
    return()
}

potentialPacks <- getFiles()

#clustering
clust <- potentialPacks %>% 
  filter(stringID == 5) %>%
  filter(Genome == unique(Genome)[4]) %>%
  getOrganismKClust(.$Genome) %>%
  print()
