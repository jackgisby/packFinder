context("Clustering Functions")

data("arabidopsisThalianaRefseq")
load("data-r/consensusSeqs.rda")
load("data-r/packClusts.rda")

consensusSeqTest <- tirClust(packClusts, arabidopsisThalianaRefseq, tirLength = 25, plot = FALSE)

test_that("no errors thrown by tirClust", {
  expect_silent(tirClust(packClusts, arabidopsisThalianaRefseq, tirLength = 25, plotSavePath = "data-raw/tirRelationships.png"))
})

test_that("consensus sequences returned are as expected", {
  expect_equal(as.character(consensusSeqTest), as.character(consensusSeqs))
})
