getPotentialPackList <- function(subSeqs, 
                                 Genome, 
                                 element.length, 
                                 TSD.length) {
  # gets potentialPack dataframe list for DNAStringSet
  
  potentialPackList <- vector("list", length = length(subSeqs))
  
  for(subSeq in 1:length(subSeqs)) {
    potentialPackList[subSeq] <- packSearch(subSeq = subSeqs[[subSeq]], 
                                            Genome, 
                                            mismatch = as.integer(subSeqs@ranges@NAMES[subSeq]), 
                                            element.length = element.length, 
                                            TSD.length = TSD.length)
  }
  
  return(potentialPackList)
}

savePotentialPacks <- function(potentialPackList, subSeqs) {
  #saves each dataframe to a csv
  
  for(i in 1:length(potentialPackList)) {
    write.csv(potentialPackList[[i]], 
              file = paste0("Data/Output/algorithmAssessment/report_", 
                            as.character(subSeqs[[i]]), 
                            "_", 
                            subSeqs@ranges@NAMES[i],
                            "mismatch"))
  }
}

saveKnownCacta <- function() {
  
}

saveOverallReport <- function(subSeqs, runTimes, mode = "normal") {
  
}
