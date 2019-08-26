library(Biostrings)

source("R/Search_Functions.R")


packSearch <- function(subSeq, Genome, mismatch = 0, element.length, TSD.length) {
  # ---input---
  # SearchString = DNAString object to be searched for
  # Genome = DNAStringSet object to search
  # mismatch = the maximum edit distance to be considered (indels + substitions)
  #   
  # ---returns---
  # locations of matches as a dataframe
  # a summary of the search, with match statistics
  
  #perform initial search for TIR matches
  forwardMatches <- as.data.frame(GRanges())
  forwardMatches <- identifyTIRMatches(forwardMatches, subSeq, Genome, mismatch = mismatch, strand = "+")
  
  reverseMatches <- as.data.frame(GRanges())
  reverseMatches <- identifyTIRMatches(reverseMatches, reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-")
  
  #determine potential transposable elements based on following elements
  potentialPack <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length, TSD.length)
  
  return(potentialPack)
}