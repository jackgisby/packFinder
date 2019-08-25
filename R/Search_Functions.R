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

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, subSeq, Genome, mismatch, element.length, TSD.length) {
  
  potTransposons <- GRanges()
  
  for(forwardMatch in 1:length(forwardMatches)) {
    
    forwardRepeat <- forwardMatches[forwardMatch]
    chr <- as.character(seqnames(forwardRepeat)@values)
    searchRange <- forwardMatches[forwardMatch]@ranges@start
    searchRange <- c(searchRange + element.length[1], searchRange + element.length[2])
    
    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }
    
    reverseRepeats <- reverseMatches[seqnames(reverseMatches) == seqnames(forwardRepeat)@values
                                     & end(reverseMatches) > searchRange[1]
                                     & end(reverseMatches) < searchRange[2]
                                     & strand(reverseMatches) == "-"]
    
    if(length(reverseRepeats) > 1) {
      
      for(reverseMatch in 1:length(reverseRepeats)) {
        fTSD <- flank(forwardRepeat, TSD.length)@ranges #todo: will cause a bug if flank is out of bounds
        rTSD <- flank(reverseRepeats[reverseMatch], TSD.length)@ranges
        
        if(Genome[Genome@ranges@NAMES == chr][[1]][fTSD] == Genome[Genome@ranges@NAMES == chr][[1]][rTSD]) {
          
          potTransposons <- c(potTransposons, GRanges(seqnames = chr,
                                                    ranges = IRanges(start = fTSD@start, 
                                                                     end = rTSD@start + rTSD@width - 1)))
        }
      }
    }
  }
  return(potTransposons)
}

potentialPacks <- identifyPotentialPackElements(forwardMatches, reverseMatches, subSeq, Genome, 1, c(300, 5000), 3)
potentialPacks
