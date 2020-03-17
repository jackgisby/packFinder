context("collapse")

data(packMatches)
data(arabidopsisThalianaRefseq)

packMatches$start <- 1
packMatches$end <- 10

collapsed <- collapseSeqs(packMatches, arabidopsisThalianaRefseq)

test_that("matches are collapsed", {
    expect_equal(collapsed$seqnames, as.factor("Chr3"))
    expect_equal(collapsed$start, 1)
    expect_equal(collapsed$end, 10)
    expect_equal(collapsed$width, 10)
    expect_equal(nrow(collapsed), 1)
})

unlink("tempTestOutput", recursive = TRUE)
