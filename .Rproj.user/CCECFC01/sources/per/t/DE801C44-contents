library(Biostrings)

source("R/Search_Functions.R")


packSearch <- function(subSeq, Genome, mismatch = 0, element.length = c(300, 5000)) {
  # ---input---
  # SearchString = DNAString object to be searched for
  # Genome = DNAStringSet object to search
  # mismatch = the maximum edit distance to be considered (indels + substitions)
  #   
  # ---returns---
  # locations of matches as a dataframe
  # a summary of the search, with match statistics
  
  #perform initial search for TIR matches
  forwardMatches <- GRanges()
  forwardMatches <- identifyTIRMatches(forwardMatches, subSeq, Genome, mismatch = mismatch, strand = "+")
  
  reverseMatches <- GRanges()
  reverseMatches <- identifyTIRMatches(reverseMatches, reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-")
  
  #determine potential transposable elements based on following elements
  potentialPack <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, element.length)
  
  return(TIR_Matches)
}