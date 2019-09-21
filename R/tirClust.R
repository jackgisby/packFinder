#' @title Analyse TIR Sequences of Pre-clustered Transposable Elements
#' @description Takes transposable elements clustered by VSEARCH, \code{\link{packClust}}, and produces consensus sequences for the terminal inverted repeats of each. Allows for the visualisation of TIR similarities between clusters for both forward and reverse strands.
#' @param packMatches Dataframe containing potential packTYPE elements and cluster information.
#' @param Genome The DNAStringSet used to create the \code{packMatches} dataframe. Required for TIR extraction.
#' @param plot Argument specifying whether the TIR consensus sequences should be plottted as a dendrogram.
#' @param plotSavePath File path for the dendrogram plot. If unspecified, the dendrogram plot is not saved.
#' @param k The k-mer size to be used for calculating a distance matrix between TIR consensus sequences. See \code{kdistance::kmer}.
#' @param tirLength The TIR size to be considered. Consensus sequences will be generated based on the first and last \code{tirLength} bases.
#' @param output Controls the output of \code{tirClust}.
#' @author Jack Gisby
#' @return If \code{output == "consensus"}, returns a list of consensus sequences for each cluster specified in \code{packMatches}. Else if \code{output == "dendrogram"}, returns a dendrogram object used to create hierarchical clustering diagrams.
#' @export

tirClust <- function(packMatches,
                     Genome,
                     plot = TRUE,
                     plotSavePath = NULL,
                     k = 5,
                     tirLength = 25,
                     output = "consensus") {
  fConsensusSeqs <- vector("list",
    length = length(unique(packMatches$clustID))
  )
  rConsensusSeqs <- vector("list",
    length = length(unique(packMatches$clustID))
  )

  for (c in 1:length(unique(packMatches$cluster))) {
    clustID <- unique(packMatches$cluster)[c]
    clust <- packMatches[packMatches$cluster == clustID, ]
    forwardTirs <- vector("list", length = length(clust[, 1]))
    reverseTirs <- vector("list", length = length(clust[, 1]))

    for (i in 1:length(clust[, 1])) {
      if (clust$strand[i] == "+") {
        forwardTirs[[i]] <- Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][clust$start[i]:(clust$start[i] + tirLength)]
        reverseTirs[[i]] <- Biostrings::reverseComplement(Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][(clust$end[i] - tirLength):clust$end[i]])
      } else if (clust$strand[i] == "-") {
        forwardTirs[[i]] <- Biostrings::reverseComplement(Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][(clust$end[i] - tirLength):clust$end[i]])
        reverseTirs[[i]] <- Genome[Genome@ranges@NAMES == clust$seqnames[i]][[1]][clust$start[i]:(clust$start[i] + tirLength)]
      }
    }

    fConsensusSeqs[[c]] <- Biostrings::consensusString((Biostrings::DNAStringSet(forwardTirs)))
    rConsensusSeqs[[c]] <- Biostrings::consensusString((Biostrings::DNAStringSet(reverseTirs)))
  }

  fConsensusSeqs <- Biostrings::DNAStringSet(unlist(fConsensusSeqs))
  rConsensusSeqs <- Biostrings::DNAStringSet(unlist(rConsensusSeqs))
  fConsensusSeqs@ranges@NAMES <- paste0("f", as.character(unique(packMatches$cluster)))
  rConsensusSeqs@ranges@NAMES <- paste0("r", as.character(unique(packMatches$cluster)))
  consensusSeqs <- c(fConsensusSeqs, rConsensusSeqs)

  dend <- stats::as.dendrogram(stats::hclust(kmer::kdistance(ape::as.DNAbin(consensusSeqs), k = k)))


  if (plot == TRUE) {
    plot(dend, main = "Clustered Transposon TIR Relationships")
  }

  if (!is.null(plotSavePath)) {
    grDevices::png(plotSavePath, width = 1000, height = 600)
    plot(dend, main = "Clustered Transposon TIR Relationships")
    grDevices::dev.off()
  }

  if (output == "dendrogram") {
    return(dend)
  } else {
    return(consensusSeqs)
  }
}
