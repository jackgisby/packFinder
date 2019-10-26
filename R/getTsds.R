#' @title
#'     Get Flanking Terminal Site Duplication Sequences
#'
#' @description
#'     Gets the flanking TSD sequences of TIRs or 
#'     predicted Pack-TYPE transposable elements. 
#'     A dataframe of these elements can be in 
#'     \code{tirMatches}.
#'
#' @param tirMatches
#'     A dataframe containing genomic ranges and 
#'     names referring to TIR sequences or predicted 
#'     Pack-TYPE transposable elements. Should be 
#'     in the format used by \code{\link{packSearch}}.
#'
#' @param Genome
#'     A DNAStringSet object containing sequences 
#'     referred to in \code{tirMatches}.
#'
#' @param tsdLength
#'     The length of the TSD region to be retrieved (integer).
#'
#' @param strand
#'     The strand of the TIR; "+" for forward, "-" 
#'     for reverse. If the TSD sequences of 
#'     transposable elements are being predicted, 
#'     then this parameter can be left as default 
#'     ("+"); if the TSD sequences of TIRs are 
#'     being found then the strand direction must 
#'     be supplied.
#'
#' @param output
#'     The type of object to be returned:
#'     \itemize{
#'         \item output = "DNAStringSet", returns a
#'             \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} 
#'             object.
#'         \item output = "character", returns a
#'             code{character} vector (default).
#'     }
#'
#' @details
#'     Called by \code{\link{packSearch}}. It is 
#'     recommended to use the general pipeline 
#'     function \code{\link{packSearch}} for 
#'     identification of potential pack elements, 
#'     which returns TSD sequences as a feature 
#'     of results, however each stage may be 
#'     called individually.
#'
#' @return
#'     Flanking TSD sequences as a vector of 
#'     characters, or if output is specified as 
#'     "DNAStringSet", TSD sequences will be 
#'     returned as a
#'     \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} 
#'     object.
#' 
#' @examples 
#'     data(arabidopsisThalianaRefseq)
#'     data(packMatches)
#' 
#'     tsdSeqs <- getTsds(packMatches, arabidopsisThalianaRefseq, 3)
#' 
#' @author
#'     Jack Gisby
#'
#' @export

getTsds <- function(tirMatches, Genome, tsdLength, strand = "+",
                        output = "character") {
    if (strand != "+" & strand != "-") {
        stop("Argument 'strand' must be specified as a character, '+' or '-'")
    }

    if (output != "character" & output != "DNAStringSet") {
        stop("Argument 'output' must be specified 
            as 'string' or 'DNAStringSet'")
    }

    if (strand == "-") {
        TSDs <- mapply(function(seqnames, end, tsdLength, Genome) {
                seq <- Genome[Genome@ranges@NAMES == seqnames][[1]]
                return(as.character(seq[(end + 1):(end + tsdLength)]))
            },
            tirMatches$seqnames,
            tirMatches$end,
            MoreArgs = list(tsdLength = tsdLength, Genome = Genome)
        )
        return(TSDs)
    } else if (strand == "+") {
        TSDs <- mapply(function(seqnames, start, tsdLength, Genome) {
                seq <- Genome[Genome@ranges@NAMES == seqnames][[1]]
                return(as.character(seq[(start - tsdLength):(start - 1)]))
            },
            tirMatches$seqnames,
            tirMatches$start,
            MoreArgs = list(tsdLength = tsdLength, Genome = Genome)
            )

        if (output == "DNAStringSet") {
            return(Biostrings::DNAStringSet(TSDs))
        } else {
            return(TSDs)
        }
    }
}
