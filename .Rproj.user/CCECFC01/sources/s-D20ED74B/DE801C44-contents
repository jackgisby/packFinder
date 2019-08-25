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
  TIR_Matches <- GRanges()
  TIR_Matches <- identifyTIRMatches(TIR_Matches, subSeq, Genome, mismatch = mismatch, strand = "+")
  TIR_Matches <- identifyTIRMatches(TIR_Matches, reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-")
  
  #determine potential transposable elements
  #potentialPack <- identifyPotentialPackElements(subSeq, Genome, element.length, TIR_Matches)
  
  return(TIR_Matches)
}