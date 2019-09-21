context("Clustering Functions")

load("tests/testthat/data-r/consensusSeqs.rda")
load("tests/testthat/data-r/packClusts.rda")

consensusSeqTest <- tirClust(packClusts, tirLength = 25, plot = FALSE)

test_that("no errors thrown by tirClust", {
  expect_silent(tirClust(packClusts, tirLength = 25, plotSavePath = "tests/testthat/data-raw/tirRelationships.png"))
})

test_that("consensus sequences returned are as expected", {
  expect_equal(as.character(consensusSeqTest), as.character(consensusSeqs))
})
