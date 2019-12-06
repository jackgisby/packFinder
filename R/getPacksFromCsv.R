#' @title
#' Retrieve Saved packFinder Results (.csv)
#'
#' @description
#' Retrieves a dataframe of potential Pack-TYPE elements, 
#' previously saved using \code{\link{packSearch}} followed 
#' by \code{\link{packsToCsv}}.
#'
#' @param file
#' File path to predicted transposons in CSV format.
#'
#' @return
#' Dataframe in the format used by \code{\link{packSearch}}.
#'
#' @seealso 
#' \code{\link{packsToCsv}}, \code{\link[utils]{read.table}}, 
#' \code{\link{packSearch}}
#'
#' @examples
#' data(packMatches)
#' 
#' packMatches <- getPacksFromCsv(
#'     system.file("extdata", "packMatches.csv", package = "packFinder")
#' )
#' 
#' @author
#' Jack Gisby
#' 
#' @export

getPacksFromCsv <- function(file) {
    checkPermissions(file)
    
    return(utils::read.csv(file, stringsAsFactors = FALSE))
}
