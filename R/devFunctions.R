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

algorithmAssessment <- function(potentialPacks, Genome) {
  # Assesses the error rate of the Pack-TYPE transposon finding algorithm
  #
  # ---input---
  # potentialPacks: a list of identified potential transposons
  # Genome: a DNAStringSet object containing the genome being searched
  #
  # ---returns---
  # prints: error rate of algorithm based on known transposons
  # returns: a list of correctly identified transposons
  
  knownCACTA <- getArAthCACTA(Genome, "complete")
  
  identifiedCACTA <- filter(potentialPacks, 
                            potentialPacks$start %in% knownCACTA$start &
                              potentialPacks$end %in% knownCACTA$end) #does not consider chromosome
  
  #number identified
  print(paste0("Correct packCACTA identified in Arabidopsis thalania: ", 
               (length(identifiedCACTA[,1])),
               "/",
               length(knownCACTA[,1])))
  #detection rate
  print(paste0("packCACTA detection rate: ", 
               round((length(identifiedCACTA[,1])/length(knownCACTA[,1])) * 100, 2),
               "%"))
  #error rate
  print(paste0("Algorithm error rate: ", 
               round((1-(length(identifiedCACTA[,1])/length(potentialPacks[,1]))) * 100, 2),
               "%"))
  
  return(identifiedCACTA)
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
  badMatches <- which(!assessSubSeq(subSeq, DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR)), mismatch))
  badCACTA <- knownCACTA[0,]
  
  for(bad in 1:length(badMatches)) {
    if(badMatches[bad] > 10) {
      badMatches[bad] <- badMatches[bad] - 10
    }
    
    badCACTA <- rbind(badCACTA, knownCACTA[badMatches[bad],])
  }
  return(badCACTA)
}

getBadSeqs <- function(knownCACTA, subSeq, mismatch) {
  badMatches <- which(!assessSubSeq(subSeq, DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR)), mismatch))
  return(knownTIRs[badMatches])
}