
context("amatch: osa")

test_that("simple test and multiple edge cases",{
  expect_equal(amatch("aa",c("ba","bb"), method="osa"), 1L)
  expect_equal(amatch("","", method="osa"), 1L)
  expect_equal(amatch(NA,"a", method="osa"), NA_integer_)
  expect_equal(amatch(NA,"a", method="osa", nomatch=0L),0L)
  expect_equal(amatch(NA,NA, method="osa"), 1L)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="osa",matchNA=FALSE, nomatch=0L), 0L)
  expect_equal(amatch("aa","bb", method="osa",maxDist=1), NA_integer_)
})


context("amatch: Damerau-Levenshtein")

test_that("simple test and multiple edge cases",{
  expect_equal(amatch("aa", c("ba","bb"), method="dl"), 1L)
  expect_equal(amatch("","", method="dl"), 1L)
  expect_equal(amatch(NA,"a", method="dl"), NA_integer_)
  expect_equal(amatch(NA,"a", method="dl",nomatch=0L), 0L)
  expect_equal(amatch(NA,NA, method="dl"), 1L)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE), NA_integer_)
  expect_equal(amatch(NA,NA, method="dl",matchNA=FALSE,nomatch=0L), 0L)
  expect_equal(amatch("aa","bb", method="dl",maxDist=1), NA_integer_)
})


