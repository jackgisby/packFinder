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
  potentialPacks <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length, TSD.length)
  
  print("Filtering complete")
  return(potentialPacks)
}