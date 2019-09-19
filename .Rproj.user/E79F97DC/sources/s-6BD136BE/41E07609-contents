identifyTIRMatches <- function(tirMatches, subSeq, Genome, mismatch, strand = "*") {
  # searches genome for potential TIRs based on sequence similarity
  #
  # ---input---
  # tirMatches: dataframe of previously identified matches
  # subSeq: DNAString object containing the sub-sequence to be searched for
  # Genome: DNAStringSet object containing the genome to be searched
  # mismatch: numeric, acceptable edit distance between subSeq and a given sequence in Genome
  # strand: string, direction of match - forward (+) or reverse (-)
  #
  # ---returns---
  # tirMatches: dataframe of previously identified matches and matches identified during this search
  
  
  for(i in 1:length(Genome)) { 
    matches <- matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE)
    
    if(length(matches) > 0) {
      tirMatches <- rbind(tirMatches, data.frame(seqnames = names(Genome)[i], 
                                                   start = matches@ranges@start, 
                                                   end = matches@ranges@start + matches@ranges@width - 1, 
                                                   width = matches@ranges@width, 
                                                   strand = strand))
    }
  }
  
  return(tirMatches)
}

getTSDs <- function(tirMatches, Genome, tsdLength, direction) {
  # gets the TSD sequences for potential TIR ends
  #
  # ---input---
  # tirMatches: dataframe of potential TIR ends
  # Genome: DNAStringSet containing the genome to be searched
  # tsdLength: numeric, length of TSD region to be matched
  # direction: string, forward (+) or reverse (-) strand
  #
  # ---returns---
  # tirMatches dataframe with an additional column for each TIR's associated TSD sequence as a string
  
  if(direction == "+") {
    return(tirMatches %>% 
             filter(start > tsdLength) %>%
             mutate(TSD = mapply(function(seqnames, start, tsdLength, Genome) {
        return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(start - tsdLength):(start - 1)]))},
        seqnames,
        start,
        MoreArgs = list(tsdLength = tsdLength, Genome = Genome))
    ))
  } else if(direction == "-") {
      return(tirMatches %>% 
               mutate(removeMatch = mapply(function(end, seqnames, tsdLength, Genome)  {
                 return((end + tsdLength) > Genome[Genome@ranges@NAMES == seqnames][[1]]@length)},
                 end,
                 seqnames,
                 MoreArgs = list(tsdLength, Genome))) %>%
               filter(removeMatch == FALSE) %>%
               select(-c(removeMatch)) %>%
               mutate(TSD = mapply(function(seqnames, end, tsdLength, Genome) {
                 return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + tsdLength)]))},
                 seqnames,
                 end,
                 MoreArgs = list(tsdLength = tsdLength, Genome = Genome))))
  }
}

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, Genome, elementLength) {
  # identifies potential Pack-TYPE transposons by identifying neighbouring forward and reverse sequences
  # and matching them together given similar neighbouring TSD sequences
  #
  # ---input---
  # forwardMatches: dataframe containing potential forward TIR sequences
  # reverseMatches: dataframe containing potential reverse TIR sequences
  # Genome: DNAStringSet object containing the genome being searched 
  # elementLength: vector of minimum and maximum length of Pack-TYPE transposons being searched for
  # 
  # ---returns---
  # dataframe of potential Pack-TYPE transposable elements
  
  packMatches <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), strand = character())
  
  for(forwardMatch in 1:length(forwardMatches[,1])) {
    
    forwardRepeat <- forwardMatches[forwardMatch,]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- c(forwardRepeat$start + elementLength[1], forwardRepeat$start + elementLength[2])
    
    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }
    
    reverseRepeats <- filter(reverseMatches,
                             seqnames == as.character(forwardRepeat$seqnames) & 
                             end > searchRange[1] & 
                             end < searchRange[2] & 
                             strand == "-"  &
                             TSD == as.character(forwardRepeat$TSD))

    if(length(reverseRepeats[,1]) > 0) {
      for(reverseMatch in 1:length(reverseRepeats[,1])) {
        packMatches <- rbind(packMatches, data.frame(seqnames = forwardRepeat$seqnames,
                                                             start = forwardRepeat$start,
                                                             end = reverseRepeats[reverseMatch,]$end,
                                                             width = reverseRepeats[reverseMatch,]$end - forwardRepeat$start,
                                                             strand = "*"))
      }
    }
  }
  return(packMatches)
}

getSeqs <- function(packMatches, Genome) {
  packMatches %>%
    mutate(seq = mapply(function(start, end, seqnames, Genome) {
      return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][start:end]))},
    start,
    end,
    seqnames,
    MoreArgs = list(Genome = Genome))) %>%
    return()
}
