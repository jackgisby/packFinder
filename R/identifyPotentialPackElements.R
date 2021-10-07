#' @title
#' Pack Element Filtering
#'
#' @description
#' Primary filtering stage for the \code{packSearch} 
#' algorithm. Identifies potential Pack-TYPE transposable 
#' elements based on proximity of matching inverted repeats 
#' and equality of TSD sequences.
#'
#' @param forwardMatches
#' A dataframe containing genomic ranges and names referring 
#' to forwards-facing TIR sequences and their respective 
#' TSD sequences.
#'
#' @param reverseMatches 
#' A dataframe containing genomic ranges and names referring 
#' to reverse-facing TIR sequences and their respective 
#' TSD sequences.
#'
#' @param Genome 
#' A DNAStringSet object containing the matches referred to 
#' in \code{forwardMatches} and \code{reverseMatches}
#'
#' @param elementLength
#' A vector of two integers containing the minimum and 
#' maximum transposable element length.
#' 
#' @param tsdMismatch
#' An integer referring to the allowable mismatch 
#' (substitutions or indels) between a transposon's TSD
#' sequences. \code{\link[Biostrings]{matchPattern}} from Biostrings 
#' is used for pattern matching. 
#'
#' @details
#' Used by \code{\link{packSearch}} as a primariy filtering 
#' stage. Identifies matches likely to be transposons based 
#' on their TIR region, from \code{\link{identifyTirMatches}}, 
#' and their TSD region, from \code{\link{getTsds}}. It is 
#' recommended to use the general pipeline function 
#' \code{\link{packSearch}} for identification of potential 
#' pack elements, however each stage may be called 
#' individually. Note that only exact TSD matches are 
#' considered, so supplying long sequences for TSD elements 
#' may lead to false-negative results.
#'
#' @return
#' A dataframe, \code{packMatches}, containing the locations 
#' of potential Pack-TYPE transposable elements in \code{Genome}.
#' 
#' @seealso 
#' \code{packSearch}
#' 
#' @examples 
#' data(arabidopsisThalianaRefseq)
#' 
#' forwardMatches <- identifyTirMatches(
#'     Biostrings::DNAString("CACTACAA"),
#'     arabidopsisThalianaRefseq,
#'     tsdLength = 3,
#'     strand = "+"
#' )
#' 
#' reverseMatches <- identifyTirMatches(
#'     Biostrings::reverseComplement(Biostrings::DNAString("CACTACAA")),
#'     arabidopsisThalianaRefseq,
#'     tsdLength = 3,
#'     strand = "-"
#' )
#' 
#' packMatches <- identifyPotentialPackElements(
#'     forwardMatches, 
#'     reverseMatches,
#'     arabidopsisThalianaRefseq, 
#'     c(300, 3500)
#' )
#' 
#' @author
#' Jack Gisby
#' 
#' @export

identifyPotentialPackElements <- function(forwardMatches, reverseMatches, 
                                            Genome, elementLength, 
                                            tsdMismatch = 0) {
    packMatches <- initialisePackMatches()

    # for each forward match, consider/filter each nearby reverse match
    for (forwardMatch in seq_len(length(forwardMatches[, 1]))) {
        forwardRepeat <- forwardMatches[forwardMatch, ]
        chr <- as.character(forwardRepeat[[1]])
        searchRange <- c(forwardRepeat$start + elementLength[1], 
                        forwardRepeat$start + elementLength[2])

        if (searchRange[2] > length(Genome[names(Genome) == chr][[1]])) {
            searchRange[2] <- length(Genome[names(Genome) == chr][[1]])
        }
    
        reverseRepeats <- filterTsdMatches(reverseMatches, forwardRepeat, 
                                            tsdMismatch, searchRange)

        if (length(reverseRepeats[, 1]) > 0) {
            
            packMatches <- rbind(
                packMatches,
                data.frame(
                    seqnames = forwardRepeat$seqnames,
                    start = forwardRepeat$start,
                    end = min(reverseRepeats$end),
                    width = min(reverseRepeats$end) - forwardRepeat$start + 1,
                    strand = "*"
                )
            )
        }
    }
    return(packMatches)
}


filterTsdMatches <- function(reverseMatches, forwardRepeat, tsdMismatch, 
                                searchRange) {
    
    # filters for exact matches (faster), or uses matchpattern
    if (tsdMismatch == 0) {
        reverseRepeats <- reverseMatches[
            reverseMatches$seqnames == as.character(forwardRepeat$seqnames) &
                reverseMatches$end > searchRange[1] &
                reverseMatches$end < searchRange[2] &
                reverseMatches$strand == "-" &
                reverseMatches$TSD == as.character(forwardRepeat$TSD),
            ]
    } else {
        reverseRepeats <- reverseMatches[
                reverseMatches$seqnames == as.character(forwardRepeat$seqnames)
                & reverseMatches$end > searchRange[1] &
                reverseMatches$end < searchRange[2] &
                reverseMatches$strand == "-",
            ]

        if (length(reverseRepeats[, 1]) > 0) {
            reverseRepeats <- reverseRepeats[
            Biostrings::vcountPattern(as.character(forwardRepeat$TSD),
                                        reverseRepeats$TSD,
                                        max.mismatch = tsdMismatch,
                                        with.indels = TRUE) > 0,
            ]
        }
    }

    return(reverseRepeats)
}
