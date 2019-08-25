identifyTIRMatches <- function(TIR_Matches, subSeq, Genome, mismatch, strand = "*") {
  
  for(i in 1:length(Genome)) {
    matchData <- GRanges(seqnames = names(Genome)[i], 
                         ranges = as(matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE), "IRanges"),
                         strand = strand)
    
    if(!is.null(matchData)) {
      TIR_Matches <- c(TIR_Matches, matchData)
    }
  }
  
  return(TIR_Matches)
}

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length) {
  
  for(forwardMatch in 1:length(forwardMatches)) {
    searchRange <- forwardMatches[forwardMatch]@ranges@start
    searchRange <- c(searchRange + element.length[1], searchRange + element.length[2])
    
    forwardRepeat <- forwardMatches[forwardMatch]
    reverseRepeats <- reverseMatches[seqnames(reverseMatches) == seqnames(forwardRepeat)@values
                                     & end(reverseMatches) > searchRange[1]
                                     & end(reverseMatches) < searchRange[2]
                                     & strand(reverseMatches) == "-"]
    print(reverseRepeats)
  }
  
  return(potTransposons)
}