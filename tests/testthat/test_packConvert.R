context("Conversion functions")

uc <- readUc("data-raw/packMatches.uc")

test_that("object type is as expected", {
  expect_type(uc, "list")
})

test_that("object dimensions are as expected", {
  expect_equal(nrow(uc), 46)
  expect_equal(ncol(uc), 8)
  expect_equal(max(uc$cluster), 16)
})
