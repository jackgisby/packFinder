#' @title
#' Sample packFinder Output
#'
#' @description
#' A sample output from \code{\link{packSearch}} with 
#' cluster information. This dataframe is in the format produced by 
#' coercing a \code{link[GenomicRanges:GRanges-class]{GRanges}} 
#' object to a dataframe: \code{data.frame(GRanges)}. 
#'
#' @format
#' A dataframe of 9 obs. and 7 variables.
#'
#' @usage
#' data(packMatches)
#'
#' @details
#' Was obtained from running \code{\link{packSearch}} 
#' on the Arabidopsis thaliana chromosome 3 reference 
#' sequence, followed by clustering using 
#' \code{\link{packClust}}. Contains the following features:
#' \itemize{
#'     \item start - the predicted element's start base 
#'     sequence position.
#'     \item end - the predicted element's end base 
#'     sequence position.
#'     \item seqnames - character string referring to the 
#'     sequence name in \code{Genome} to which \code{start} 
#'     and \code{end} refer to.
#' }
#' 
#' The dataset was generated as in the example below.
#'
#' @examples
#' data(arabidopsisThalianaRefseq)
#' 
#' packMatches <- packSearch(
#'     Biostrings::DNAString("CACTACAA"),
#'     arabidopsisThalianaRefseq,
#'     elementLength = c(300, 3500),
#'     tsdLength = 3
#' )
#' 
#' @seealso 
#' \code{\link{packSearch}}, \code{\link{data.frame}}, 
#' \code{\link{arabidopsisThalianaRefseq}}

"packMatches"
