#' @title
#'     Retrieve Saved packFinder Results (.csv)
#'
#' @description
#'     Retrieves a dataframe of potential Pack-TYPE
#'     elements, previously saved using 
#'     \code{\link{packSearch}} followed by 
#'     \code{\link{packsToCsv}}.
#'
#' @param file
#'     path to predicted transposons in CSV format.
#'
#' @author
#'     Jack Gisby
#'
#' @return
#'     Dataframe in the format used by \code{\link{packSearch}}.
#'
#' @seealso 
#'     \code{\link{packsToCsv}}
#'
#' @examples
#'     data(packMatches)
#' 
#'     packsToCsv(packMatches, "packMatches.csv")
#'     packMatches <- getPacksFromCsv("packMatches.csv")
#' 
#' @export

getPacksFromCsv <- function(file) {
    if (!is.null(file)) {
        if (!(file.access(file, 4) == 0) |
            !(file.access(file, 4) == 0) |
            !(file.access(file, 2) == 0)) {
            stop("file does not exist, or R does not 
                have read/write permissions")
        }
    }

    return(utils::read.csv(file, stringsAsFactors = FALSE))
}
