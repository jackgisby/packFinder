library(Biostrings)

source("R/searchFunctions.R")

packSearch <- function(subSeq, Genome, mismatch = 0, element.length, TSD.length) {
  # General use pipeline function for the Pack-TYPE transposon finding algorithm
  #
  # ---input---
  # subSeq: DNAString object containing the TIR sequence
  # Genome: DNAStringSet object to be searched
  # mismatch: the maximum edit distance to be considered for TIR matches (indels + substitions)
  # element.length: vector containing the minimum and maximum length of the Pack-TYPE transposons
  #   
  # ---returns---
  # a dataframe of potential Pack-TYPE transposons
  
  #perform initial search for TIR matches and get related TSD sequences
  print("Getting forward matches")
  forwardMatches <- as.data.frame(GRanges()) %>%
    identifyTIRMatches(subSeq, Genome, mismatch = mismatch, strand = "+") %>%
    getTSDs(Genome, TSD.length, direction = "+")
  
  print("Getting reverse matches")
  reverseMatches <- as.data.frame(GRanges()) %>%
    identifyTIRMatches(reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-") %>%
    getTSDs(Genome, TSD.length, direction = "-")

  #determine potential transposable elements based on nearby elements and TSD sequences
  print("Filtering matches based on TSD sequences")
  potentialPacks <- identifyPotentialPackElements(forwardMatches, reverseMatches, Genome, element.length)
  
  print("Filtering complete")
  return(potentialPacks)
}

blastFilter <- function(potentialPacks) {
  
}