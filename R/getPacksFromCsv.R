getPacksFromCsv <- function(file) {

  if (!is.null(file)) {
    if (!(file.access(file, 4) == 0) |
        !(file.access(file, 4) == 0) |
        !(file.access(file, 2) == 0)) {
      stop("file does not exist, or R does not have read/write permissions")
    }
  }
  return(utils::read.csv(file, stringsAsFactors = FALSE))
}
