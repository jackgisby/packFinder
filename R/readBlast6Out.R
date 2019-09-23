readBlast6Out <- function(file) {

  if (!is.null(file)) {
    if (!(file.access(file, 4) == 0) |
        !(file.access(file, 4) == 0) |
        !(file.access(file, 2) == 0)) {
      stop("file does not exist, or R does not have read/write permissions")
    }
  }

  blast6out <- utils::read.table(file, sep = "\t")
  colnames(blast6out) <- c(
    "query",
    "subject",
    "identity",
    "width",
    "mismatch",
    "openGaps",
    "queryStart",
    "queryEnd",
    "targetStart",
    "targetEnd",
    "e",
    "bits"
  )

  return(blast6out)
}

