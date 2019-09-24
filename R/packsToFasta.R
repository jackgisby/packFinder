#' @title
#' Save packFinder Results in FASTA Format (.fasta)
#'
#' @description
#' Saves a dataframe of potential Pack-TYPE elements, usually generated via
#' \code{\link{packSearch}}. May be retrieved using
#' \code{\link{getPacksFromFasta}}.
#'
#' @param file
#' FASTA file save path.
#'
#' @param packMatches
#' A dataframe containing genomic ranges and names referring
#' to sequences to be extracted. Can be obtained from \code{\link{packSearch}}
#' or generated from a \code{\link[GenomicRanges]{GRanges}} object, after
#' conversion to a dataframe. Must contain the following features:
#' \itemize {
#'   \item start - the predicted element's start base sequence position.
#'   \item end - the predicted element's end base sequence position.
#'   \item seqnames - character string referring to the sequence name in
#'   \code{Genome} to which \code{start} and \code{end} refer to.
#' }
#'
#' @param Genome
#' A DNAStringSet object containing sequences referred to in \code{packMatches}
#' (the object originally used to predict the transposons
#' \code{\link{packSearch}}).
#'
#' @examples \dontrun {
#' data(arabidopsisThalianaRefseqSubset)
#'
#' packMatches <- packsToFasta(
#'   packMatches,
#'   "path/to/packMatches.csv",
#'   arabidopsisThalianaRefseq)
#' }
#'
#' @author
#' Jack Gisby
#'
#' @seealso \code{\link{getPacksFromFasta}}
#'
#' @export

packsToFasta <- function(packMatches, file, Genome) {
  file.create(file)
  for (match in 1:length(packMatches[, 1])) {
    write(c(
      paste0(
        ">", packMatches$seqnames[match],
        " | start = ", packMatches$start[match],
        " | end = ", packMatches$end[match],
        " | width = ", packMatches$width[match],
        " | strand = ", packMatches$strand[match],
        " | TSD = ", packMatches$TSD[match]
      ),
      as.character(Genome[Genome@ranges@NAMES == packMatches$seqnames[match]][[1]][packMatches$start[match]:packMatches$end[match]])
    ),
    file = file,
    append = TRUE
    )
  }
  return(print(paste0("FASTA written to ", file)))
}
