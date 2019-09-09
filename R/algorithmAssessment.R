source("R/packSearch.R")
source("R/devFunctions.R")

#functions to assess algorithm/filtering performance

assessPotentialPackList <- function(subSeqs, Genome, element.length, TSD.length) {
  # assesses each item of a potentialPackList
  
  potentialPackList <- getPotentialPackList(subSeqs, Genome, element.length, TSD.length)
}

assessRepeatMapFilter <- function() {
  # assesses a potentialPackList using a repeat map filtering stage
  
  repeatMaps <- getRepeatMaps(Genome)
  potentialPacks <- filterElements(potentialPacks, repeatMaps)
  
  knownCACTA <- saveReport(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
  knownCACTA <- saveReport(filter(potentialPacks, isTransposon == FALSE), subSeq, Genome, integrityFilter = NULL, mismatch = max.mismatch)
  print(end-start)
}

assessBlastFilter <- function() {
  # assesses a potentialPackList using a blast filtering stage
  
  #Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
  #                        "C:\\Users\\jackg\\Documents\\R\\nt_db\\ncbi-blast-2.9.0+\\bin"
  #                        , sep= .Platform$path.sep)) #sets path of blast+ exe files
  
  db <- blast(db="C:/Users/jackg/Documents/R/nt_db/nt/nt", type = "blastn")
  

}