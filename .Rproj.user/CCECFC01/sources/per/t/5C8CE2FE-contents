algorithmAssessment <- function(transposonList, Genome) {
  knownCACTA <- read.csv("Input/knownCACTA.csv", sep = ";")
  knownCACTA <- knownCACTA[,1:8]
  x <- data.frame(name = Genome@ranges@NAMES)
  
  for(i in 1:length(knownCACTA$Chr)) {
    knownCACTA$chrName[i] <- as.character(x$name[knownCACTA$Chr[i]])
  }
  
  knownCACTA$Identified <- knownCACTA$start %in% transposonList$start
  
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", sum(knownCACTA$Identified)))
  print(paste0("Algorithm error rate: ", (1-(sum(knownCACTA$Identified)/length(transposonList[,1])))))
}

algorithmAssessment(potentialPacks, Genome)
