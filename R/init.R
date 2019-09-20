#' @title Initialise Global Variables
#' @description Initialise global variables in order to prevent R variable binding issues (e.g. for \code{dplyr} & \code{mapply})
#' @note Called upon package build only.
#' @author Jack Gisby

init <- function() {
  utils::globalVariables(c(
    "Genome",
    "TSD",
    "cluster",
    "desc",
    "end",
    "getTSDs",
    "query",
    "removeMatch",
    "seqnames",
    "start",
    "strand",
    "type",
    "width"
  ))
}

init()
