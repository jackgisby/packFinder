#' @title
#' Export packFinder Results to a GRanges Object
#'
#' @description
#' A dataframe containing genomic ranges and names referring 
#' to sequences to be extracted, likely obtained from 
#' \code{\link{packSearch}}, can be converted to a GRanges 
#' object. Can be converted back to a dataframe using 
#' \code{\link{getPacksFromGRanges}}. Additional features, 
#' such as clusters and TSD sequences, will be included in 
#' the object as metadata columns.
#'
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names 
#' referring to sequences to be extracted. Can be obtained 
#' from \code{\link{packSearch}} or generated from a 
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} object, 
#' after conversion to a dataframe. Must contain the 
#' following features:
#' \itemize{
#'     \item start - the predicted element's start 
#'     base sequence position.
#'     \item end - the predicted element's end base 
#'     sequence position.
#'     \item seqnames - character string referring 
#'     to the sequence name in \code{Genome} to which 
#'     \code{start} and \code{end} refer to.
#' }
#'
#' @return
#' A GRanges object containing the ranges contained in 
#' \code{packMatches} and additional metadata columns. May 
#' be easily converted between dataframe and GRanges format 
#' for use in the \code{packFinder} package and 
#' \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' package. Note that most functions in the \code{packFinder} 
#' package require sequence ranges to be provided in 
#' dataframe format.
#'
#' @seealso 
#' \code{\link{getPacksFromGRanges}}
#'     
#' @author
#' Jack Gisby
#' 
#' @examples 
#' data(packMatches)
#' packGRanges <- packsToGRanges(packMatches)
#'
#' @export

packsToGRanges <- function(packMatches) {
    return(GenomicRanges::GRanges(packMatches))
}
