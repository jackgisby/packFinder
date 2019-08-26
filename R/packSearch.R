library(Biostrings)

source("R/Search_Functions.R")

#overall pipeline for pack finding algorithm
packSearch <- function(subSeq, Genome, mismatch = 0, element.length, TSD.length) {
  # ---input---
  # SearchString = DNAString object to be searched for
  # Genome = DNAStringSet object to search
  # mismatch = the maximum edit distance to be considered (indels + substitions)
  #   
  # ---returns---
  # locations of matches as a dataframe
  # a summary of the search, with match statistics
  
  #perform initial search for TIR matches and get related TSD sequences
  forwardMatches <- as.data.frame(GRanges()) %>%
    identifyTIRMatches(subSeq, Genome, mismatch = mismatch, strand = "+") %>%
    getTSDs(Genome, direction = "+")
  
  reverseMatches <- as.data.frame(GRanges()) %>%
    identifyTIRMatches(reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-") %>%
    getTSDs(Genome, direction = "-")
    
  
  #determine potential transposable elements based on nearby elements and TSD sequences
  potentialPack <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length, TSD.length)
  
  return(forwardMatches)#potentialPack
}