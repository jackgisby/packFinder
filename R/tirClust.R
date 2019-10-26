#' @title Analyse TIR Sequences of Pre-clustered Transposable Elements
#'
#' @description
#' Takes transposable elements clustered by VSEARCH, \code{\link{packClust}},
#' and produces consensus sequences for the terminal inverted repeats of each.
#' Allows for the visualisation of TIR similarities between clusters for both
#' forward and reverse strands.
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names referring
#' to sequences to be extracted. Can be obtained from \code{\link{packSearch}}
#' or generated from a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object, after
#' conversion to a dataframe. Must contain the following features:
#' \itemize{
#'   \item start - the predicted element's start base sequence position.
#'   \item end - the predicted element's end base sequence position.
#'   \item seqnames - character string referring to the sequence name in
#'   \code{Genome} to which \code{start} and \code{end} refer to.
#' }
#'
#' @param Genome
#' A DNAStringSet object containing sequences referred to in \code{packMatches}
#' (the object originally used to predict the transposons
#' \code{\link{packSearch}}).
#'
#' @param plot
#' Argument specifying whether the TIR consensus sequences should be plottted
#' as a dendrogram.
#'
#' @param plotSavePath
#' File path for the dendrogram plot. If unspecified, the dendrogram plot is
#' not saved.
#'
#' @param k
#' The k-mer size to be used for calculating a distance matrix between TIR
#' consensus sequences. See \code{kdistance::kmer}. Larger word sizes will not
#' be suitable for longer TIR sequences, due to processing time required.
#' Additionally k must be greater than the TIR sequence length.
#'
#' @param tirLength
#' The TIR size to be considered. Consensus sequences will be generated based
#' on the first and last \code{tirLength} bases of a transposon.
#'
#' @param output
#' Controls the output of \code{tirClust}. If output is specified as
#' "consensus", the consensus sequences of each TIR cluster will be returned;
#' else, if output is specified as "dendrogram", a dendrogram object will be
#' returned for creation of customisable plots.
#'
#' @author
#' Jack Gisby
#'
#' @return
#' If \code{output} is specified as "consensus" (default), returns a list of
#' consensus sequences for each cluster specified in \code{packMatches} as a
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}. Else if \code{output} is specified
#' as "dendrogram", returns a dendrogram object used to create hierarchical
#' clustering diagrams.
#'
#' @seealso
#' code{\link{packClust}}, code{\link{packAlign}}
#'
#' @export

tirClust <- function(packMatches, Genome, tirLength = 25, plot = TRUE,
                     plotSavePath = NULL, k = 5, output = "consensus") {
    if (output != "consensus" & output != "dendrogram") {
        stop("Argument 'output' must be specified as 'consensus' or 'dendrogram'")
    }

    fConsensusSeqs <- vector("list", length = length(unique(packMatches$clustID)))
    rConsensusSeqs <- vector("list", length = length(unique(packMatches$clustID)))

    for (c in seq_len(length(unique(packMatches$cluster)))) {
        clustID <- unique(packMatches$cluster)[c]
        clust <- packMatches[packMatches$cluster == clustID, ]
        forwardTirs <- vector("list", length = length(clust[, 1]))
        reverseTirs <- vector("list", length = length(clust[, 1]))

        for (i in seq_len(length(clust[, 1]))) {
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
        grDevices::png(plotSavePath, width = 1500, height = 1000)
        plot(dend, main = "Clustered Transposon TIR Relationships")
        grDevices::dev.off()
    }

    if (output == "dendrogram") {
        return(dend)
    } else {
        return(consensusSeqs)
    }
}
