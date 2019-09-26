context("data")

data("arabidopsisThalianaRefseq")
data("packMatches")

test_that("length of datasets is as expected", {
  expect_equal(length(arabidopsisThalianaRefseq), 1)
  expect_equal(sum(arabidopsisThalianaRefseq@ranges@width), 3800001)
  expect_equal(length(packMatches[,1]), 6)
  expect_equal(length(packMatches[1,]), 7)
})

test_that("type of dataset is as expected", {
  expect_is(arabidopsisThalianaRefseq, "DNAStringSet")
  expect_is(packMatches, "data.frame")
})
