getPotentialPackList <- function(subSeqs, 
                                 Genome, 
                                 element.length, 
                                 TSD.length) {
  # gets potentialPack dataframe list for DNAStringSet
  
  potentialPackList <- NULL
  runTimes <- vector("list", length = length(subSeqs))
  
  for(subSeq in 1:length(subSeqs)) {
    start = Sys.time()
    potentialPacks <- packSearch(subSeq = subSeqs[[subSeq]], 
                                            Genome = Genome, 
                                            mismatch = as.integer(subSeqs@ranges@NAMES[subSeq]), 
                                            element.length = element.length, 
                                            TSD.length = TSD.length)
                                        
    potentialPacks$stringID <- subSeq
    runTimes[[subSeq]] <- Sys.time() - start
    
    if(!is.null(potentialPackList)) {
      potentialPackList <- rbind(potentialPackList, potentialPacks)
    } else {
      potentialPackList <- potentialPacks
    }
  }
  
  return(list(potentialPackList = potentialPackList,
              runTimes = runTimes))
}

saveOverallReport <- function(subSeqs, runTimes, mode, detectRate = NULL, errorRate = NULL, errorTotal) {
  # saves overall report
  
  if(mode == "Arath") {
    overallReport <- data.frame(Search_ID = 1:length(subSeqs),
                                Search_Sequence = as.character(subSeqs),
                                Allowable_Mismatch = subSeqs@ranges@NAMES,
                                Run_Time = unlist(runTimes),
                                Detected = detectRate,
                                Errors = errorRate,
                                Total_Matches = errorTotal)
  } else {
    overallReport <- data.frame(Search_ID = 1:length(subSeqs),
                                Search_Sequence = as.character(subSeqs),
                                Allowable_Mismatch = subSeqs@ranges@NAMES,
                                Run_Time = unlist(runTimes),
                                Total_Matches = errorTotal)
  }
  
  write.csv(overallReport, "Data/Output/algorithmAssessment/overallReport.csv", row.names = FALSE)
}

saveKnownCacta <- function(subSeqs, potentialPackList, Genome, runTimes, integrityFilter = NULL) {
  # saves known cacta as report and makes specialised overall report
  
  knownCactas <- NULL
  detectRate <- vector("integer", length(subSeqs))
  errorRate <- vector("integer", length(subSeqs))
  errorTotal <- vector("integer", length(subSeqs))
  
  for(subSeq in 1:length(subSeqs)) {
    potentialPacks <- filter(potentialPackList, stringID == subSeq)
    knownCACTA <- getArathCACTA(Genome, integrityFilter)  %>%
      mutate(forwardTIR_Identified = assessSubSeq(subSeqs[[subSeq]], getknownTIRs(.), mismatch = as.integer(subSeqs@ranges@NAMES[subSeq]))[1:10]) %>%
      mutate(reverseTIR_Identified = assessSubSeq(subSeqs[[subSeq]], getknownTIRs(.), mismatch = as.integer(subSeqs@ranges@NAMES[subSeq]))[11:20]) %>%
      mutate(forwardTIR = mapply(function(forwardTIR) {
        return(as.character(forwardTIR))},
        forwardTIR)) %>% 
      mutate(reverseTIR = mapply(function(reverseTIR) {
        return(as.character(reverseTIR))},
        reverseTIR)) %>%
      mutate(Chr = mapply(function(Chr) {
        return(Genome@ranges@NAMES[Chr])},
        Chr)) %>%
      select(-c(TAIR10.annotations, blast.best.hits, mobilization, chrNames)) %>%
      mutate(Search_ID = subSeq) %>%
      mutate(identified = start %in% potentialPacks$start 
             & end %in% potentialPacks$end
             & Chr %in% potentialPacks$seqnames)
    
    if(!is.null(knownCactas)) {
      knownCactas <- rbind(knownCactas, knownCACTA)
    } else {
      knownCactas <- knownCACTA
    }
    
    detectRate[subSeq] <-sum(knownCACTA$identified)
    errorRate[subSeq] <- length(potentialPacks[,1]) - sum(knownCACTA$identified)
    errorTotal[subSeq] <- length(potentialPacks[,1])
  }
  
  write.csv(knownCactas, "Data/Output/algorithmAssessment/knownCACTA.csv", row.names = FALSE)
  
  saveOverallReport(subSeqs, runTimes, mode = "Arath", detectRate, errorRate, errorTotal)
}

getArathCACTA <- function(Genome, integrityFilter = NULL) {
  # gets the ArAth packCACTA sequences
  #
  # ---input---
  # integrityFilter: (optional) string, filters knownCACTA - "complete" filters for 
  # only complete matches whereas "not partial" filters for non-partial matches
  # Genome: DNAStringSet object containing the ArAth genome
  #
  # ---returns---
  # dataframe containing sequence information from the known ArAth CACTA sequences
  
  knownCACTA <- read.csv("Data/Data/knownCACTA.csv", sep = ";") %>%
    mutate(TSD = gsub("\\*", "", TSD)) %>%
    mutate(chrNames = Genome@ranges@NAMES[Chr]) %>%
    mutate(forwardTIR = mapply(function(Chr, start, Genome) {
      return(Genome[[Chr]][start:(start+25)])},
      Chr,
      start,
      MoreArgs = list(Genome = Genome))) %>%
    mutate(reverseTIR = mapply(function(Chr, end, Genome) {
      return(reverseComplement(Genome[[Chr]][(end-25):end]))},
      Chr,
      end,
      MoreArgs = list(Genome = Genome)))
  
  if (is.null(integrityFilter)) {
    return(knownCACTA)
    
  } else if (integrityFilter == "complete") {
    filter(knownCACTA, integrity == integrityFilter) %>%
      return()
    
  } else if (integrityFilter == "not partial") {
    filter(knownCACTA, integrity != "partial") %>%
      return()
  }
}

getknownTIRs <- function(knownCACTA) {
  # gets the TIRs of the known CACTA sequences as a DNAStringSet
  
  return(DNAStringSet(c(knownCACTA$forwardTIR, knownCACTA$reverseTIR)))
}

assessSubSeq <- function(subSeq, knownTIRs, mismatch = 0) {
  # assesses a given subSeq for recognition of the known TIR sequences
  
  successfulMatches <- vector(mode = "logical", length = length(knownTIRs))
  
  for(i in 1:length(knownTIRs)) {
    if(countPattern(subSeq, knownTIRs[[i]], max.mismatch = mismatch, with.indels = TRUE) > 0) {
      successfulMatches[i] <- TRUE
    }
  }
  
  return(successfulMatches)
}
