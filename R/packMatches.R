#' @title Sample packFinder Output
#' @format A dataframe containing 9 obs. of 7 variables.
#' @docType data
#' @usage data(packMatches)
#' @details The sample data was generated as in the example, using the Arabidopsis thaliana chromosome 3 subset data.
#' @examples
#' \dontrun{
#' data(arabidopsisThalianaRefseq)
#'
#' subSeq <- Biostrings::DNAString("CACTACAA")
#' }
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
#'   identity = 0.5
#' )
#'

"packMatches"
