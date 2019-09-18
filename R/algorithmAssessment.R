#pipeline functions to assess algorithm/filtering performance
assessPotentialPackList <- function(subSeqs, Genome, element.length, TSD.length, mode = "normal") {
  # assesses each item of a potentialPackList
  
  potentialPackListRuntimes <- getPotentialPackList(subSeqs = subSeqs, 
                                            Genome = Genome, 
                                            element.length = element.length, 
                                            TSD.length = TSD.length)
  potentialPackList <- potentialPackListRuntimes[[1]]
  runTimes <- potentialPackListRuntimes[[2]]
  unlink("Data/Output/algorithmAssessment/", recursive = TRUE)
  write.csv(potentialPackList, "Data/Output/algorithmAssessment/potentialPacks.csv", row.names = FALSE)
  print("potentialPack report saved")
  
  if(mode == "Arath") {
    saveKnownCacta(subSeqs, potentialPackList, Genome, runTimes)
    print("knownCACTA report saved")
    print("Overall report saved")
  } else {
    errorTotal <- vector("integer", length(subSeqs))
    for(i in 1:max(potentialPackList$stringID)) {
      errorTotal[i] <- filter(potentialPackList, stringID == i) %>%
        select(stringID) %>%
        unlist() %>%
        length()
    }
    saveOverallReport(subSeqs, runTimes, mode = "normal", errorTotal = errorTotal)
    print("Overall report saved")
  }
}


assessRepeatMapFilter <- function() {
  # assesses a potentialPackList using a repeat map filtering stage
  
  
}

assessBlastFilter <- function() {
  # assesses a potentialPackList using a blast filtering stage
  
  #Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
  #                        "C:\\Users\\jackg\\Documents\\R\\nt_db\\ncbi-blast-2.9.0+\\bin"
  #                        , sep= .Platform$path.sep)) #sets path of blast+ exe files
  
  db <- blast(db="C:/Users/jackg/Documents/R/nt_db/nt/nt", type = "blastn")
  

}
