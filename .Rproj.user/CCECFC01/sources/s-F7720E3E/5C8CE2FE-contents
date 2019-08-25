algorithmAssessment <- function(transposonList) {
  knownCACTA <- read.csv("Input/knownCACTA.csv", sep = ";")
  score <- 0
  
  for(element in 1:length(knownCACTA)) {
    i = which(transposonList$start == knownCACTA$start[element])
    if(length(i) == 0) {
      i <- 0
    } else { #if (transposonList$end[i] == knownCACTA$end[element]) {
      score <- score + 1
    }
  }
  
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", score))
  print(paste0("Algorithm error rate: ", (1-(score/length(transposonList[,1])))))
}

algorithmAssessment(potentialPacks)
