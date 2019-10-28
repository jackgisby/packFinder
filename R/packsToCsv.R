#' @title
#' Save packFinder Results in CSV Format (.csv)
#'
#' @description
#' Saves a dataframe of potential Pack-TYPE elements, 
#' usually generated via \code{\link{packSearch}}. May be 
#' retrieved using \code{\link{getPacksFromCsv}}.
#'
#' @param file
#' CSV file save path.
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names referring 
#' to sequences to be extracted. Can be obtained from 
#' \code{\link{packSearch}} or generated from a 
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} object, 
#' after conversion to a dataframe. Must contain the 
#' following features:
#' \itemize{
#'     \item start - the predicted element's start base 
#'     sequence position.
#'     \item end - the predicted element's end base 
#'     sequence position.
#'     \item seqnames - character string referring to 
#'     the sequence name in \code{Genome} to which 
#'     \code{start} and \code{end} refer to.
#' }
#'
#' @return
#' Save location of csv file.
#'
#' @seealso 
#' \code{\link{getPacksFromCsv}}, \code{\link[utils]{write.table}},
#' \code{\link{packSearch}}
#' 
#' @examples
#' data(packMatches)
#' 
#' packsToCsv(
#'     packMatches, 
#'     system.file("extdata", "packMatches.csv", package = "packFinder")
#' )
#' 
#' @author
#' Jack Gisby
#'
#' @export

packsToCsv <- function(packMatches, file) {
    utils::write.csv(packMatches, file, row.names = FALSE)

    return(print(paste0("File successfully written to ", file)))
}
