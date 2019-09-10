# useful functions for manipulating and assessing other functions

getGenomeDnaStringSet <- function(genomeName = "Arabidopsis thaliana", genomePath = "Input/") {
  # Loads the ArAth genome and required packages for testing
  #
  # ---returns---
  # Arabidopsis thalania genome (as Biostrings::DNAStringSet)
  
  Genome <- read_genome(getGenome(db = "refseq", genomeName, path = genomePath))
  if(genomeName == "Arabidopsis thaliana") {
    return(Genome[1:5])
  } else {
    return(Genome)
  }
}

getRepeatMaps <- function(Genome) {
  # gets map of Arath repeats
  
  repeatMaps <- data.frame("Chromosome" = factor(),
                           "Start" = integer(),
                           "End" = integer(), 
                           "Name" = character(), 
                           "Rep_Start" = integer(), 
                           "Rep_End" = integer(), 
                           "Orientation" = factor(), 
                           "Identity_(%)" = double())
  
  chr <- c("I", "II", "III", "IV", "V")
  
  for(i in 1:length(Genome@ranges@NAMES)) {
    repeatMap <- read.table(paste0("Input/AraTh_RepeatMap/ATmap", chr[[i]]))
    colnames(repeatMap) <- c("Chromosome", "Start", "End", "Name", "Rep_Start", "Rep_End", "Orientation", "Identity_(%)")
    repeatMap$Chromosome <- Genome@ranges@NAMES[i]
    repeatMaps <- rbind(repeatMaps, repeatMap)
  }
  
  repeatMaps %>%
    filter(grepl("ENSPM", Name)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    return()
}

getFastaFromDataFrame <- function(dataframe, Genome, filepath) {
  # saves dataframe as fasta file
  
  getDNAStringSetFromDataFrame(dataframe, Genome) %>%
    writeXStringSet(filepath)
  
  print(paste0("FASTA file successfully written to: ", filepath))
}
