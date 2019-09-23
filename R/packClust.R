#' @title Cluster Transposons using VSEARCH
#' @description Cluster potential pack-TYPE elements by sequence similarity.
#' @param packMatches A dataframe of potential transposable elements. Will be saved as a FASTA file for VSEARCH.
#' @param Genome The \code{\link[Biostrings]{DNAStringSet}} object used to create the \code{packMatches} dataframe. Used to retrieve sequences for clustering.
#' @param identity The sequence identity of two transposable elements in \code{packMatches} required to be grouped into a cluster.
#' @param threads The number of threads to be used by VSEARCH.
#' @param strand The strand direction (+, - or *) to be clustered.
#' @param saveFolder The folder to save output files (uc, blast6out, FASTA)
#' @param vSearchPath The location of the VSEARCH executable file.
#' @note
#' In order to cluster sequences using VSEARCH, the executable file must first be installed.
#' @author Jack Gisby
#' @return Saves cluster information, including a \code{uc} and \code{blast6out} file, to the specified location. Returns the given \code{packMatches} dataframe with an additional column, \code{cluster}, containing cluster IDs.
#' @references VSEARCH may be downloaded from \url{https://github.com/torognes/vsearch}. See \url{https://www.ncbi.nlm.nih.gov/pubmed/27781170} for further information.
#' @export



packClust <- function(packMatches,
                      Genome,
                      identity = 0.6,
                      threads = 1,
                      strand = "both",
                      saveFolder = NULL,
                      vSearchPath = "path/to/vsearch/vsearch-2.14.1-win-x86_64/vsearch.exe") {
  if (is.null(saveFolder)) {
    saveFolder <- getwd()
  }

  if (parallel::detectCores() < threads) {
    stop("There are not ", threads, " cores available")
  }

  if (identity > 1 | identity < 0) {
    stop("Identity must be of type integer or double, and have a value between 0 and 1")
  }

  if (strand == "+") {
    strand <- "plus"
  } else if (strand == "*") {
    strand <- "both"
  } else if (strand != "both" & strand != "plus") {
    message("Argument 'strand' must be specified as 'plus' or 'both'")
  }

  if (system2(vSearchPath, "--v") != 0) {
    message("VSEARCH cannot be found at", vSearchPath, ". Please ensure that
            the VSEARCH executable file is installed at the correct location
            and that the full file path is specified.")
  }

  packMatchesFile <- paste0(saveFolder, "packMatches.fasta")
  ID <- as.integer(rownames(packMatches))
  packMatches$ID <- ID
  packMatches <- packMatches[order(-packMatches$width), ]

  packMatchesSet <- getPackSeqs(packMatches, Genome, output = "DNAStringSet")
  packMatchesSet@ranges@NAMES <- as.character(rownames(packMatches))
  Biostrings::writeXStringSet(packMatchesSet, packMatchesFile)

  system2(
    command = vSearchPath,
    args = paste0(
      "--cluster_smallmem ", packMatchesFile, " \ ",
      "--id ", identity, " \ ",
      "--strand ", strand, " \ ",
      "--threads ", threads, " \ ",
      "--qmask none \ ",
      "--uc ", file.path(saveFolder, paste0("packMatches", ".uc")), " \ ",
      "--log ", file.path(saveFolder, paste0("packMatches", ".log")), " \ ",
      "--blast6out ", file.path(saveFolder, paste0("packMatches", ".blast6out")), " \ ",
      "--clusterout_sort \ ",
      "--clusterout_id \ ",
      "--sizeout"
    )
  )

  vSearchClusts <- readUc(file = file.path(saveFolder, paste0("packMatches", ".uc")))
  vSearchClusts <- vSearchClusts[vSearchClusts$type != "C", ]

  packMatches$strand <- mapply(function(strand) {
    if (strand == "*") {
      return("+")
    } else {
      return(strand)
    }
  },
  strand = as.character(vSearchClusts$strand)
  )

  packMatches$cluster <- vSearchClusts$cluster
  rownames(packMatches) <- packMatches$ID
  packMatches <- packMatches[order(packMatches$ID), ]
  packMatches <- subset(packMatches, select = -c(ID))
  return(packMatches)
}
