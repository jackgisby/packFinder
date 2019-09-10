#pipeline functions to assess algorithm/filtering performance

assessPotentialPackList <- function(subSeqs, Genome, element.length, TSD.length, mode = "normal") {
  # assesses each item of a potentialPackList
  
  potentialPackList <- getPotentialPackList(subSeqs, Genome, element.length, TSD.length)
  
  unlink("Data/Output/algorithmAssessment/", recursive = TRUE)
  
  savePotentialPacks(potentialPackList)
  
  if(mode == "Arath") {
    saveKnownCacta()
    saveOverallReport(mode = "Arath")
  } else {
    saveOverallReport(mode = "normal")
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