context("Conversion functions")

ucClust <- readUc("data-raw/clustOutputuc")
ucAlign <- readUc("data-raw/alignOutput.uc")
blast6Out <- readBlast6Out("data-raw/clustOutput.blast6out")

test_that("object type is as expected", {
  expect_type(ucClust, "list")
  expect_type(ucAlign, "list")
  expect_type(blast6Out, "list")
})

test_that("object dimensions are as expected", {
  expect_equal(nrow(ucClust), 46)
  expect_equal(ncol(ucClust), 8)

  expect_equal(nrow(uc), 46)
  expect_equal(ncol(uc), 6)

  expect_equal(nrow(blast6Out), 46)
  expect_equal(ncol(blast6Out), 8)
})

test_that("object output is as expected", {
  expect_equal(max(ucClust$cluster), 16)
})
