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
  forwardMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character()) %>%
    identifyTIRMatches(subSeq, Genome, mismatch = mismatch, strand = "+") %>%
    getTSDs(Genome, TSD.length, direction = "+")
  
  print("Getting reverse matches")
  reverseMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character()) %>%
    identifyTIRMatches(reverseComplement(subSeq), Genome, mismatch = mismatch, strand = "-") %>%
    getTSDs(Genome, TSD.length, direction = "-")
  
  #case: no matches
  if(length(forwardMatches[,1]) == 0 | length(reverseMatches[,1]) == 0) {
    return(None)
  }
  
  #determine potential transposable elements based on nearby elements and TSD sequences
  print("Filtering matches based on TSD sequences")
  potentialPacks <- identifyPotentialPackElements(forwardMatches, reverseMatches, Genome, element.length) %>%
    getTSDs(Genome, TSD.length, "+") %>%
    getTIRs(Genome)
  
  print("Initial filtering complete")
  return(potentialPacks)
}

packBlast <- function(potentialPacks, db, db.loc = "local", Genome) {
  blastMatches <- getBlastMatches(potentialPacks, db, db.loc) %>%
    return()
}

packPipeline <- function(subSeq, Genome, mismatch = 0, element.length, TSD.length, db.loc = "online") {
  packSearch(subSeq, Genome, mismatch, element.length, TSD.length) %>%
    packBlast() %>%
    return()
}
