packSearch <- function(subSeq, Genome, mismatch = 0, elementLength, tsdLength) {
  # General use pipeline function for the Pack-TYPE transposon finding algorithm
  #
  # ---input---
  # subSeq: DNAString object containing the TIR sequence
  # Genome: DNAStringSet object to be searched
  # mismatch: the maximum edit distance to be considered for TIR matches (indels + substitions)
  # elementLength: vector containing the minimum and maximum length of the Pack-TYPE transposons
  #   
  # ---returns---
  # a dataframe of potential Pack-TYPE transposons
  
  #perform initial search for TIR matches and get related TSD sequences
  print("Getting forward matches")
  forwardMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character()) %>%
    identifyTIRMatches(subSeq, Genome, mismatch = mismatch, strand = "+") %>%
    getTSDs(Genome, tsdLength, direction = "+")
  
  print("Getting reverse matches")
  reverseMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character()) %>%
    identifyTIRMatches(reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-") %>%
    getTSDs(Genome, tsdLength, direction = "-")
  
  #case: no matches
  if(length(forwardMatches[,1]) == 0 | length(reverseMatches[,1]) == 0) {
    return(None)
  }
  
  #determine potential transposable elements based on nearby elements and TSD sequences
  print("Filtering matches based on TSD sequences")
  packMatches <- identifyPotentialPackElements(forwardMatches, reverseMatches, Genome, elementLength) %>%
    getTSDs(Genome, tsdLength, "+") %>%
    getSeqs(Genome)
  
  print("Initial filtering complete")
  return(packMatches)
}
