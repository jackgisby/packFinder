#' @title
#' Remove Low Quality Sequences
#'
#' @description
#' Takes transposable elements detected by 
#' \code{\link{packSearch}} and removes those with large 
#' numbers of wildcard ("N") bases. Used by 
#' \code{\link{packClust}} and \code{\link{packAlign}} to 
#' remove poor quality sequences that may interfere with 
#' the quality of sequence alignments.
#'
#' @param packMatches 
#' A dataframe containing genomic ranges 
#' and names referring to sequences to be extracted.
#'
#' @param Genome
#' The original set of sequences used to generate the 
#' transposons detected by \code{\link{packSearch}}.
#'
#' @param maxWildcards
#' The maximal allowable proportion of wildcards in the 
#' sequence of each match (defaults to \code{0.05}).
#' 
#' @return
#' The original dataframe, \code{packMatches}, with 
#' sequences removed that are found to contain a proportion 
#' of wildcards ("N") greater than that specified in 
#' \code{maxWildcards}.
#' 
#' @author
#' Jack Gisby
#' 
#' @examples
#' data(arabidopsisThalianaRefseq)
#' data(packMatches)
#'     
#' filteredMatches <- filterWildcards(
#'     packMatches, 
#'     arabidopsisThalianaRefseq, 
#'     maxWildcards = 0.05
#' )
#' 
#' @seealso 
#' \code{\link{packClust}}, \code{\link{packAlign}}, \code{packMatches},
#' data(arabidopsisThalianaRefseq)
#' 
#' @author 
#' Jack Gisby
#' 
#' @export

filterWildcards <- function(packMatches, Genome, maxWildcards = 0.05) {
    badMatches <- vector(mode = "logical", length = nrow(packMatches))
    
    message(paste0("Removing matches with a proportion of 
                    wildcards ('N's) above ", maxWildcards))
    
    # gets each sequence and removes those with many "N" bases
    for (i in seq_len(nrow(packMatches))) {
        seq <- Genome[names(Genome) == packMatches$seqnames[i]]
        seq <- as.character(seq[[1]][packMatches$start[i]:packMatches$end[i]])
        
        if (grepl("N", seq)) {
            if ((nchar(gsub("N", "", seq)) / nchar(seq)) > maxWildcards) {
                badMatches[i] <- TRUE
            }
        }
    }

    message(paste0(sum(badMatches), " matches removed"))
    return(packMatches[badMatches == FALSE, ])
}
