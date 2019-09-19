# useful functions for manipulating and assessing other functions

getGenomeDnaStringSet <- function(genomeName = "Arabidopsis thaliana", 
                                  genomePath = "Data/Data/genomes/",
                                  db = "refseq") {
  # Loads the ArAth genome and required packages for testing
  #
  # ---returns---
  # Arabidopsis thalania genome (as Biostrings::DNAStringSet)
  
  Genome <- read_genome(getGenome(db = db, genomeName, path = genomePath, reference = TRUE))
  # if(genomeName == "Arabidopsis thaliana") {
  #   return(Genome[1:5])
  # } else {
  #   return(Genome)
  # }
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


getFiles <- function() {
  #get previously generated data from folders
  files <- vector("list", length(list.files("Data/Output/algorithmAssessment/")))
  
  for(file in 1:length(list.files("Data/Output/algorithmAssessment/"))) {
    files[[file]] <- read.csv(paste0("Data/Output/algorithmAssessment/",
                                     list.files("Data/Output/algorithmAssessment/")[file],
                                     "/potentialPacks.csv"))
    
    files[[file]]$Genome <- list.files("Data/Output/algorithmAssessment/")[file]
    
  }
  
  rbind(files[[1]],
        files[[2]],
        files[[3]],
        files[[4]]) %>%
    filter(TSD != "NNN") %>%
    return()
}

knownCacta <- function(Genome) {
  knownCACTA <- read.csv("Data/Data/knownCACTA.csv", stringsAsFactors = FALSE, sep = ";")
  
  data.frame(start = knownCACTA$start,
             end = knownCACTA$end,
             width = knownCACTA$end - knownCACTA$start,
             strand = "*") %>%
    mutate(seqnames = Genome@ranges@NAMES[as.integer(knownCACTA$Chr)]) %>%
    getTSDs(Genome, direction = "+", tsdLength = 3) %>%
    getSeqs(Genome) %>%
    return()
}