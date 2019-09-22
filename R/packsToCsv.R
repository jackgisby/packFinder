packsToCsv <- function(packMatches, file) {
  write.csv(packMatches, file)

  return(print(paste0("File successfully written to ", file)))
}
