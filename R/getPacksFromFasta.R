#' @title
#' Retrieve Saved packFinder Results (.fasta)
#'
#' @description
#' Retrieves a dataframe of potential Pack-TYPE elements, 
#' previously saved using \code{\link{packSearch}} followed 
#' by \code{\link{packsToFasta}}. Parses the .fasta file 
#' and title field containing:
#' \itemize{
#'     \item seqnames - name of origin sequence
#'     \item start - transposon base start position on 
#'     origin sequence
#'     \item end - transposon base end position on origin 
#'     sequence
#'     \item width - width of transposon
#'     \item strand - direction of transposon 
#'     ("+", "-" or "*")
#'     \item TSD - terminal site duplication (TSD) sequence
#' }
#'
#' @param file
#' Path to predicted transposons in FASTA format.
#'
#' @return
#' Dataframe in the format used by \code{\link{packSearch}}.
#'
#' @seealso
#' \code{\link{packsToFasta}}, \code{\link{packSearch}}
#' 
#' @examples
#' data(arabidopsisThalianaRefseq)
#' data(packMatches)
#' 
#' packMatches <- getPacksFromFasta(
#'     system.file("extdata", "packMatches.fasta", package = "packFinder")
#' )
#' 
#' @author
#' Jack Gisby
#'
#' @export

getPacksFromFasta <- function(file) {
    checkPermissions(file)
    fileCon <- file(file, "r")
    
    packMatches <- initialisePackMatches(TSD = TRUE)
    
    while (TRUE) {
        seqName <- readLines(fileCon, n = 1)
        
        # read FASTA until end of file
        if ((length(seqName) == 0) | (length(seq) == 0)) {
            break
        } else if (substr(seqName, 1, 1) != ">") {
            break
        }
        seq <- readLines(fileCon, n = 1)

        # read package-specific FASTA file into dataframe
        seqName <- gsub(">", "", seqName)
        seqName <- gsub("start =", "", seqName)
        seqName <- gsub("end =", "", seqName)
        seqName <- gsub("width =", "", seqName)
        seqName <- gsub("strand =", "", seqName)
        seqName <- gsub("TSD =", "", seqName)
        seqName <- gsub(" ", "", seqName)
        seqName <- strsplit(seqName, "|", fixed = TRUE)[[1]]
        names(seqName) <- c("seqnames", "start", "end", 
                            "width", "strand", "TSD")
        
        packMatches <- rbind(packMatches, seqName, stringsAsFactors = FALSE)
        seqName <- NULL
    }
    
    close(fileCon)
    colnames(packMatches) <- c("seqnames", "start", "end", 
                                "width", "strand", "TSD")
    return(packMatches)
}

checkPermissions <- function(file) {
    if (!is.null(file)) {
        # check for read/write permissions
        if (!(file.access(file, 2) == 0)) {
            stop("File does not exist, or R does not 
                have read permissions")
        }
    }
}

initialisePackMatches <- function(TSD = FALSE) {
    if (TSD == TRUE) {
        packMatches <- data.frame(
            seqnames = character(),
            start = integer(),
            end = integer(),
            width = integer(),
            strand = character(),
            TSD = character()
        )
    } else {
        packMatches <- data.frame(
            seqnames = character(),
            start = integer(),
            end = integer(),
            width = integer(),
            strand = character()
        )
    }
    
    return(packMatches)
}