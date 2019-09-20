read.uc <- function(savePath) {
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
