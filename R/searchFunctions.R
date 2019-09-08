identifyTIRMatches <- function(TIR_Matches, subSeq, Genome, mismatch, strand = "*") {
  # searches genome for potential TIRs based on sequence similarity
  #
  # ---input---
  # TIR_Matches: dataframe of previously identified matches
  # subSeq: DNAString object containing the sub-sequence to be searched for
  # Genome: DNAStringSet object containing the genome to be searched
  # mismatch: numeric, acceptable edit distance between subSeq and a given sequence in Genome
  # strand: string, direction of match - forward (+) or reverse (-)
  #
  # ---returns---
  # TIR_Matches: dataframe of previously identified matches and matches identified during this search
  
  
  for(i in 1:length(Genome)) { #refactor without Granges?
    matchData <- GRanges(seqnames = names(Genome)[i], 
                         ranges = as(matchPattern(subSeq, Genome[[i]], max.mismatch = mismatch, with.indels = TRUE), "IRanges"),
                         strand = strand)
    
    if(!is.null(matchData)) {
      TIR_Matches <- rbind(TIR_Matches, as.data.frame(matchData))
    }
  }
  
  return(TIR_Matches)
}

getTSDs <- function(TIR_Matches, Genome, TSD.length, direction) {
  # gets the TSD sequences for potential TIR ends
  #
  # ---input---
  # TIR_Matches: dataframe of potential TIR ends
  # Genome: DNAStringSet containing the genome to be searched
  # TSD.length: numeric, length of TSD region to be matched
  # direction: string, forward (+) or reverse (-) strand
  #
  # ---returns---
  # TIR_Matches dataframe with an additional column for each TIR's associated TSD sequence as a string
  
  if(direction == "+") {
    return(TIR_Matches %>% mutate(
      TSD = mapply(function(seqnames, start, TSD.length, Genome) {
        return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(start - TSD.length):(start - 1)]))},
        seqnames,
        start,
        MoreArgs = list(TSD.length = TSD.length, Genome = Genome))
    ))
  } else if(direction == "-") {
    return(TIR_Matches %>% mutate(
      TSD = mapply(function(seqnames, end, TSD.length, Genome) {
        return(as.character(Genome[Genome@ranges@NAMES == seqnames][[1]][(end + 1):(end + TSD.length)]))},
        seqnames,
        end,
        MoreArgs = list(TSD.length = TSD.length, Genome = Genome))))
  }
}

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, Genome, element.length) {
  # identifies potential Pack-TYPE transposons by identifying neighbouring forward and reverse sequences
  # and matching them together given similar neighbouring TSD sequences
  #
  # ---input---
  # forwardMatches: dataframe containing potential forward TIR sequences
  # reverseMatches: dataframe containing potential reverse TIR sequences
  # Genome: DNAStringSet object containing the genome being searched 
  # element.length: vector of minimum and maximum length of Pack-TYPE transposons being searched for
  # 
  # ---returns---
  # dataframe of potential Pack-TYPE transposable elements
  
  potTransposons <- as.data.frame(GRanges())
  
  for(forwardMatch in 1:length(forwardMatches[,1])) {
    
    forwardRepeat <- forwardMatches[forwardMatch,]
    chr <- as.character(forwardRepeat[[1]])
    searchRange <- c(forwardRepeat$start + element.length[1], forwardRepeat$start + element.length[2])
    
    if(searchRange[2] > length(Genome[Genome@ranges@NAMES == chr][[1]])) {
      searchRange[2] <- length(Genome[Genome@ranges@NAMES == chr][[1]])
    }
    
    reverseRepeats <- filter(reverseMatches,
                             seqnames == forwardRepeat$seqnames & 
                             end > searchRange[1] & 
                             end < searchRange[2] & 
                             strand == "-"  &
                             TSD == forwardRepeat$TSD)

    if(length(reverseRepeats[,1]) > 0) {
      for(reverseMatch in 1:length(reverseRepeats[,1])) {
        potTransposons <- rbind(potTransposons, data.frame(seqnames = forwardRepeat$seqnames,
                                                             start = forwardRepeat$start,
                                                             end = reverseRepeats[reverseMatch,]$end,
                                                             width = reverseRepeats[reverseMatch,]$end - forwardRepeat$start,
                                                             strand = "*"))
      }
    }
  }
  return(potTransposons)
}

getBlastMatches <- function(potentialPacks, db, db.loc, Genome) {
  
  
  if(db.loc == "local") {
    blastMatches <- list(length = length(potentialPacks))
    for(i in 1:length(potentialPacks)) {
      blastMatches[i] <- blastMatches[i] <- getBlastMatchesLocal(Genome[Genome@ranges@NAMES == potentialPacks$seqnames[i]][potentialPacks$start[i]:potentialPacks$end[i]],db)
    }
  } else if(db.loc == "online") {
    blastMatches <- getBlastMatchesOnline(potentialPacks, db, Genome)
  }
  return(blastMatches)
}

getBlastMatchesLocal <- function(DNAStringSetQuery, db) {
  predict(db, 
          DNAStringSetQuery,
          BLAST_args = "-num_threads 4") %>%
    return()
}

getDNAStringSetFromDataFrame <- function(dataframe, Genome) {
  strings <- DNAStringSet()
  for(i in 1:length(dataframe[,1])) {
    strings <- c(strings, DNAStringSet(Genome[Genome@ranges@NAMES == dataframe$seqnames[i]][[1]][dataframe$start[i]:dataframe$end[i]]))
  }
  
  return(DNAStringSet(strings))
}

getBlastMatchesOnline <- function(potentialPacks, db = "nt", Genome) {
  packStrings <- getDNAStringSetFromDataFrame(potentialPacks, Genome)
  
  return(blastSeq(packStrings,
                  n_blast = length(packStrings),
                  delay_req = 20,
                  database="nt",
                  keepInMemory = TRUE,
                  email = "jackgisby@gmail.com"))
}
