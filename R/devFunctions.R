# useful functions for manipulating and assessing other functions within the pack-TYPE
# transposon finding project

initialise <- function() {
  # Loads the ArAth genome and required packages for testing
  #
  # ---returns---
  # Arabidopsis thalania genome (as Biostrings::DNAStringSet)
  
  library(Biostrings)
  library(biomartr)
  library(GenomicRanges)
  library(dplyr)
  
  Genome <- read_genome(getGenome(db = "refseq", "Arabidopsis thaliana", path = "/Input"))
  return(Genome[1:5])
}

getArAthCACTA <- function(Genome, integrityFilter = NULL) {
  # gets the ArAth packCACTA sequences 
  #
  # ---input---
  # integrityFilter: (optional) string, filters knownCACTA - "complete" filters for 
  # only complete matches whereas "not partial" filters for non-partial matches
  # Genome: DNAStringSet object containing the ArAth genome
  #
  # ---returns---
  # dataframe containing sequence information from the known ArAth CACTA sequences


  
  knownCACTA <- read.csv("Input/knownCACTA.csv", sep = ";") %>%
    mutate(TSD = gsub("\\*", "", TSD)) %>%
    mutate(chrNames = Genome@ranges@NAMES[Chr]) %>%
    mutate(forwardTIR = mapply(function(Chr, start, Genome) {
      return(Genome[[Chr]][start:(start+25)])},
      Chr,
      start,
      MoreArgs = list(Genome = Genome))) %>%
    mutate(reverseTIR = mapply(function(Chr, end, Genome) {
      return(reverseComplement(Genome[[Chr]][(end-25):end]))},
      Chr,
      end,
      MoreArgs = list(Genome = Genome)))
    
    if (is.null(integrityFilter)) {
      return(knownCACTA)
      
    } else if (integrityFilter == "complete") {
      filter(knownCACTA, integrity == integrityFilter) %>%
        return()
      
    } else if (integrityFilter == "not partial") {
      filter(knownCACTA, integrity != "partial") %>%
        return()
    }
}

algorithmAssessment <- function(potentialPacks, Genome, integrityFilter = "complete") {
  # Assesses the error rate of the Pack-TYPE transposon finding algorithm
  #
  # ---input---
  # potentialPacks: a list of identified potential transposons
  # Genome: a DNAStringSet object containing the genome being searched
  # integrityFilter: (optional) string, filters knownCACTA - "complete" filters for 
  # only complete matches whereas "not partial" filters for non-partial matches
  # Genome: DNAStringSet object containing the ArAth genome
  #
  # ---returns---
  # prints: error rate of algorithm based on known transposons
  # returns: a list of correctly identified transposons
  
  knownCACTA <- getArAthCACTA(Genome, integrityFilter) %>%
    mutate(identified = start %in% potentialPacks$start & end %in% potentialPacks$end) #does not consider chromosome

  
  #number identified
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", 
               (sum(knownCACTA$identified)),
               "/",
               length(knownCACTA[,1])))
  #detection rate
  print(paste0("packCACTA detection rate: ", 
               round((sum(knownCACTA$identified)/length(knownCACTA[,1])) * 100, 2),
               "%"))
  #error rate
  print(paste0("Algorithm error rate: ", 
               round((1-(knownCACTA$identified/length(potentialPacks[,1]))) * 100, 2),
               "%"))
  
  return(knownCACTA)
}

getknownTIRs <- function(knownCACTA) {
  return(DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR)))
}

assessSubSeq <- function(subSeq, knownTIRs, mismatch = 0) {
  successfulMatches <- vector(mode = "logical", length = length(knownTIRs))
  
  for(i in 1:length(knownTIRs)) {
    if(countPattern(subSeq, knownTIRs[[i]], max.mismatch = mismatch, with.indels = TRUE) > 0) {
      successfulMatches[i] <- TRUE
    }
  }
  
  return(successfulMatches)
}

getBadMatches <- function(knownCACTA, subSeq, mismatch) {
  badMatches <- which(!assessSubSeq(subSeq, getknownTIRs(knownCACTA), mismatch))
  badCACTA <- knownCACTA[0,]
  
  for(bad in 1:length(badMatches)) {
    if(badMatches[bad] > 10) {
      badMatches[bad] <- badMatches[bad] - 10
    }
    
    badCACTA <- rbind(badCACTA, knownCACTA[badMatches[bad],])
  }
  return(badCACTA)
}



getBadSeqs <- function(subSeq, mismatch) {
  badMatches <- which(!assessSubSeq(subSeq, getknownTIRs(knownCACTA), mismatch))
  return(knownTIRs[badMatches])
}

saveReport <- function(potentialPacks, subSeq, Genome, integrityFilter = NULL, mismatch = 0) {
  knownCACTA <- algorithmAssessment(potentialPacks, Genome, integrityFilter = integrityFilter) %>%
    mutate(forwardTIR_Identified = assessSubSeq(subSeq, getknownTIRs(.), mismatch = mismatch)[1:10]) %>%
    mutate(reverseTIR_Identified = assessSubSeq(subSeq, getknownTIRs(.), mismatch = mismatch)[11:20]) %>%
    mutate(forwardTIR = mapply(function(forwardTIR) {
      return(as.character(forwardTIR))},
      forwardTIR)) %>% 
    mutate(reverseTIR = mapply(function(reverseTIR) {
      return(as.character(reverseTIR))},
      reverseTIR)) %>%
    rename(chr = Chr) %>%
    select(-c(TAIR10.annotations, blast.best.hits, mobilization, chrNames))
    
  write.csv(knownCACTA, file = "Results/algorithmReport.csv")
  return(knownCACTA)
}
