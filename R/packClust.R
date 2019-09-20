#' Clustering of potential transposable elements using VSEARCH.
#' @param packMatches A dataframe of potential transposable elements. Will be saved as a FASTA file for VSEARCH.
#' @param identity The sequence identity of two transposable elements in \code{packMatches} required to be grouped into a cluster.
#' @param threads The number of threads to be used by VSEARCH.
#' @param strand The strand direction (+, - or *) to be clustered.
#' @param saveFolder The folder to save output files (uc, blast6out, FASTA)
#' @param vSearchPath The location of the VSEARCH executable file.
#' @return Saves cluster information, including a \code{uc} and \code{blast6out} file, to the specified location. Returns the given \code{packMatches} dataframe with an additional column, \code{cluster}, containing cluster IDs.
#' @export

packClust <- function(packMatches,
                      identity = 0.6,
                      threads = 1,
                      strand = "both",
                      saveFolder = "packFinder/vSearch/",
                      vSearchPath = "path/to/vsearch/vsearch-2.14.1-win-x86_64/vsearch.exe") {


  packMatchesFile <- saveFolder
  packMatches <- packMatches %>%
    mutate(ID = 1:length(packMatches[,1])) %>%
    arrange(desc(width))

  packMatchesSet <- DNAStringSet(packMatches$seq)
  packMatchesSet@ranges@NAMES <- as.character(packMatches$ID)
  writeXStringSet(packMatchesSet, packMatchesFile)

  system2(
    command = vSearchPath,
    args = paste0(
      "--cluster_smallmem ",
      packMatchesFile,
      " \ ",
      "--qmask none \ ",
      "--uc ",
      file.path(saveFolder, paste0("packMatches", ".uc")),
      " \ ",
      "--id ",
      identity,
      " \ ",
      "--threads ",
      threads,
      " \ ",
      "--clusterout_sort \ ",
      "--clusterout_id \ ",
      "--strand ",
      strand,
      " \ ",
      "--log ",
      file.path(saveFolder, paste0("packMatches", ".log")),
      " \ ",
      "--blast6out ",
      file.path(saveFolder, paste0("packMatches", ".blast6out")),
      " \ ",
      "--sizeout"
    )
  )

  vSearchClusts <- read.uc(file.path(saveFolder, paste0("packMatches", ".uc"))) %>%
    filter(type != "C") %>%
    arrange(query)
  packMatches %>%
    mutate(strand = mapply(function(strand) {
      if(strand == "*") {
        return("+")
      } else {
        return(strand)
      }
    },
    strand = as.character(vSearchClusts$strand))) %>%
    mutate(cluster = vSearchClusts$cluster) %>%
    return()
}
