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
