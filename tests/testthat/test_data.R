context("data")

data("arabidopsisThalianaRefseq")

test_that("length of datasets is as expected", {
  expect_equal(length(arabidopsisThalianaRefseq), 7)
  expect_equal(sum(arabidopsisThalianaRefseq@ranges@width), 119668634)
})

test_that("type of dataset is as expected", {
  expect_is(arabidopsisThalianaRefseq, "DNAStringSet")
})

test_that("name of dataset is as expected", {
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[1],
               "NC_003070.9 Arabidopsis thaliana chromosome 1 sequence")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[2],
               "NC_003071.7 Arabidopsis thaliana chromosome 2 sequence")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[3],
               "NC_003074.8 Arabidopsis thaliana chromosome 3 sequence")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[4],
               "NC_003075.7 Arabidopsis thaliana chromosome 4 sequence")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[5],
               "NC_003076.8 Arabidopsis thaliana chromosome 5 sequence")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[6],
               "NC_037304.1 Arabidopsis thaliana ecotype Col-0 mitochondrion, complete genome")
  expect_equal(arabidopsisThalianaRefseq@ranges@NAMES[7],
               "NC_000932.1 Arabidopsis thaliana chloroplast, complete genome")
})
