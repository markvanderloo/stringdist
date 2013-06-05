
context("qgrams")

test_that("qgram edge cases",{
  expect_equal(qgrams('a' ,      q=1), as.table(c(a=1)))   # basic test
  expect_equal(qgrams('aa',      q=1), as.table(c(a=2)))   # idem
  expect_equal(qgrams(c('a','a'),q=1), as.table(c(a=2)))   # count unique n-grams
  expect_equal(qgrams(c(NA,'a'), q=1), as.table(c(a=1)))   # skip NA's
  expect_equal(qgrams(NA,q=1), table(integer(0)))          # skip all
  expect_equivalent(qgrams(c("a","ab"), q=2), table("ab")) # skip q>nchar
  expect_equal(qgrams(c("a"),q=2), table(integer(0)))      # skip all
  expect_equivalent(qgrams(c(''),q=0), table(''))          # empty string, q=0
})


