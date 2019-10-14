context("Clustering Functions")

dir.create("tempTestOutput")
data("arabidopsisThalianaRefseq")
data("packMatches")
load("data-r/consensusSeqs.rda")

consensusSeqTest <- tirClust(packMatches,
                             arabidopsisThalianaRefseq,
                             tirLength = 25,
                             plot = FALSE)

test_that("no errors thrown by tirClust", {
  expect_silent(tirClust(packMatches,
                         arabidopsisThalianaRefseq,
                         tirLength = 25,
                         plot = FALSE,
                         plotSavePath = "tempTestOutput/tirRelationships.png"))
})

test_that("consensus sequences returned are as expected", {
  expect_equal(as.character(consensusSeqTest), as.character(consensusSeqs))
})

unlink("tempTestOutput", recursive = TRUE)
