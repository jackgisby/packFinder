#' @title
#' Retrieve Saved packFinder Results (.fasta)
#'
#' @description
#' Retrieves a dataframe of potential Pack-TYPE elements, previously saved using
#' \code{\link{packSearch}} followed by \code{\link{packsToFasta}}.
#' Parses the .fasta file and title field containing:
#' \itemize{
#'  \item seqnames - name of origin sequence
#'  \item start - transposon base start position on origin sequence
#'  \item end - transposon base end position on origin sequence
#'  \item width - width of transposon
#'  \item strand - direction of transposon ("+", "-" or "*")
#'  \item TSD - terminal site duplication (TSD) sequence
#' }
#'
#' @param file
#' Path to predicted transposons in FASTA format.
#'
#' @examples
#' \dontrun{
#' packMatches <- getPacksFromFasta("path/to/packMatches.fasta")
#' }
#'
#' @author
#' Jack Gisby
#'
#' @return
#' Dataframe in the format used by \code{\link{packSearch}}.
#'
#' @seealso
#' \code{\link{packsToFasta}}
#'
#' @export

getPacksFromFasta <- function(file) {
    if (!is.null(file)) {
        if (!(file.access(file, 4) == 0) |
            !(file.access(file, 4) == 0) |
            !(file.access(file, 2) == 0)) {
            stop("file does not exist, or R does not 
                have read/write permissions")
        }
    }

    fileCon <- file(file, "r")
    packMatches <- data.frame(
        seqnames = character(),
        start = integer(),
        end = integer(),
        width = integer(),
        strand = character(),
        TSD = character()
    )
    while (TRUE) {
        seqName <- readLines(fileCon, n = 1)
        if ((length(seqName) == 0) | (length(seq) == 0)) {
            break
        } else if (substr(seqName, 1, 1) != ">") {
            break
        }
        seq <- readLines(fileCon, n = 1)

        seqName <- gsub(">", "", seqName)
        seqName <- gsub("start =", "", seqName)
        seqName <- gsub("end =", "", seqName)
        seqName <- gsub("width =", "", seqName)
        seqName <- gsub("strand =", "", seqName)
        seqName <- gsub("TSD =", "", seqName)
        seqName <- gsub(" ", "", seqName)
        seqName <- strsplit(seqName, "|", fixed = TRUE)[[1]]
        names(seqName) <- c("seqnames", "start", "end", "width", "strand", "TSD")
        packMatches <- rbind(packMatches,
            seqName,
            stringsAsFactors = FALSE
        )
        seqName <- NULL
    }
    close(fileCon)
    colnames(packMatches) <- c("seqnames", "start", "end", 
                                "width", "strand", "TSD")
    return(packMatches)
}
