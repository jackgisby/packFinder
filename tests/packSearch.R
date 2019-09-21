context("packSearch")
library(testthat)
library(packFinder)

Genome <- data("arabidopsisThalianaRefseq")
subSeq <- Biostrings::DNAString("CACTACAA")

packMatches <- packSearch(subSeq, Genome, mismatch = 0, elementLength = c(300, 3500), tsdLength = 3)

test_that("dimensions of packMatches are correct", {
  expect_equal(nrow(packMatches), 29)
  expect_equal(ncol(packMatches), 7)
})

