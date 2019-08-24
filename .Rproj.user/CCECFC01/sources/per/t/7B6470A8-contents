library(Biostrings)

source("R/Visualise_Chromomap.R")
source("R/Search_Biostrings.R")



Search_Visualise <- function(SearchString, Genome, max.mismatch = 0, with.indels = FALSE, fixed = TRUE) {
  # ---input---
  # SearchString = DNAString object to be searched for
  # Genome = BSgenome object to search
  # searchChr = chromosomes to search (list of strings)
  # max.mismatch = the maximum edit distance to be considered
  # with.indels = whether or not to consider indels as a mismatch
  # fixed = allows the use of IUPAC symbols
  #   
  # ---returns---
  # saves visualisation to "Plots" folder
  # saves locations of matches to tsv files
  # a summary of the search, with match statistics
  
  chrList = seqnames(Genome)
  matches = getMatches(SearchString, Genome, chrList, max.mismatch = max.mismatch, with.indels = with.indels, fixed = fixed)
  
  if(matches == "No matches found.") {
    return(matches)
  }
  
  return(Visualise(Genome, matches[[1]], matches[[2]], chrList))
}