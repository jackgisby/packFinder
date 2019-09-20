#' Function for the conversion of uc files to a dataframe.
#' @param savePath The file path of the .uc file.
#' @return A dataframe, \code{packClusts}, containing the converted .uc file.
#' @export

readUc <- function(savePath) {
  packClusts <- read.table(savePath, sep = "\t")
  colnames(packClusts) <- c("type",
                            "cluster",
                            "width",
                            "identity",
                            "strand",
                            "6",
                            "7",
                            "cigarAlignment",
                            "query",
                            "target"
  )

  packClusts %>%
    select(-c("6", "7")) %>%
    return()
}
