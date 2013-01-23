
library(testthat)

context("General ")
test_that("Argument parsing",{
   expect_equal(stringdist(character(0),"a"),numeric(0))
   expect_equal(stringdist("a",character(0)),numeric(0))
}) 


context("Optimal String Alignment")
test_that("Edge cases in OSA method",{
   expect_equal(stringdist( "","a",method='osa'),1)
   expect_equal(stringdist("a", "",method='osa'),1)
   expect_equal(stringdist("a","a",method='osa'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='osa',maxDist=1),-1)
})

test_that("weights are handled correctly",{
   # deletion
   expect_equal(stringdist("ab","a", method='osa',weight=c(0.5,1,1,1)),0.5)
   # insertion
   expect_equal(stringdist("a" ,"ab",method='osa',weight=c(1,0.5,1,1)),0.5)
   # substitution
   expect_equal(stringdist("a" ,"b", method='osa',weight=c(1,1,0.5,1)),0.5)
   # transposition
   expect_equal(stringdist("ac","ca",method='osa',weight=c(1,1,1,0.5)),0.5)
   # symmetry property in simple case
   expect_equal(
      stringdist("abc","ac",method='dl',weight=c(0.5,1,1,1)),
      stringdist("ac","abc",method='dl',weight=c(1,0.5,1,1))
   )
})

context("Damerau-Levenstein")




