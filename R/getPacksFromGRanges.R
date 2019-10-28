#' @title
#' Retrieve packFinder Results from GRanges Object
#'
#' @description
#' A \code{link[GenomicRanges:GRanges-class]{GRanges}}  
#' object, potentially 
#' generated using \code{\link{packSearch}} and 
#' \code{\link{packsToGRanges}}, can be converted to a 
#' dataframe. If a GRanges object is supplied without TSD 
#' information, this can be calculated and appended to the 
#' final dataframe.
#'
#' @param packGRanges
#' \code{link[GenomicRanges:GRanges-class]{GRanges}}
#' object to be coerced.
#'
#' @param Genome
#' (optional) Sequences referred to by \code{packGRanges}.
#'
#' @param tsdLength
#' (optional) Length of TSD sequences.
#'
#' @return 
#' Dataframe in the format used by \code{\link{packSearch}}. 
#' If \code{Genome} and \code{tsdLength} are supplied, then 
#' TSD sequences are retrieved and returned 
#' as part of the dataframe.
#' 
#' @seealso
#' \code{\link{packsToGRanges}}, 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}},
#' \code{\link{packSearch}} 
#' 
#' @examples
#' data(packMatches)
#' 
#' GRangesObject <- packsToGRanges(packMatches)
#' packMatches <- getPacksFromGRanges(GRangesObject)
#' 
#' @author
#' Jack Gisby
#'
#' @export

getPacksFromGRanges <- function(packGRanges, Genome = NULL, tsdLength = NULL) {
    if (is.null(Genome) | is.null(tsdLength)) {
        return(as.data.frame(packGRanges))
    }
    else if (!is.null(Genome) & !is.null(tsdLength)) {
        packMatches <- as.data.frame(packGRanges)
        packMatches$TSDs <- getTsds(packMatches, Genome, tsdLength)
        return(packMatches)
    }
}
