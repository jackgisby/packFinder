#' @title
#' Sample packFinder Output
#'
#' @description
#' A sample output from \code{\link{packSearch}} with cluster information.
#'
#' @format
#' A dataframe of 9 obs. and 7 variables.
#'
#' @usage
#' data(packMatches)
#'
#' @details
#' Was obtained from running \code{\link{packSearch}} on the Arabidopsis
#' thaliana chromosome 3 reference sequence, followed by clustering using
#' \code{\link{packClust}}. Contains the following features:
#' \itemize {
#'   \item start - the predicted element's start base sequence position.
#'   \item end - the predicted element's end base sequence position.
#'   \item seqnames - character string referring to the sequence name in
#'   \code{Genome} to which \code{start} and \code{end} refer to.}
#' The dataset was generated as in the example below.
#'
#' @examples
#' \dontrun{
#' data(arabidopsisThalianaRefseq)
#'
#' subSeq <- Biostrings::DNAString("CACTACAA")
#'
#' packMatches <- packSearch(subSeq,
#'   arabidopsisThalianaRefseq,
#'   mismatch = 0,
#'   elementLength = c(300, 3500),
#'   tsdLength = 3
#' )
#'
#' packMatches <- packClust(packMatches,
#'   arabidopsisThalianaRefseq,
#'   saveFolder = "devData/",
#'   vSearchPath = "D:/vsearch-2.14.1-win-x86_64/vsearch.exe",
#'   identity = 0.5,
#'   identityDefinition = 2)
#'   }

"packMatches"
