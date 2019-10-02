context("packSearch")

data("arabidopsisThalianaRefseq")
data("packMatches")

packMatchesTest <- packSearch(
  Biostrings::DNAString("CACTACAA"),
  arabidopsisThalianaRefseq,
  mismatch = 0,
  elementLength = c(300, 3500),
  tsdLength = 3
)

test_that("dimensions of packMatches are correct", {
  expect_equal(nrow(packMatchesTest), 6)
  expect_equal(ncol(packMatchesTest), 6)

  expect_equal(nrow(packMatches), 6)
  expect_equal(ncol(packMatches), 7)
})

test_that("packSearch returns correct results", {
  expect_equal(packMatchesTest, subset(packMatches, select = -c(cluster)))
})
