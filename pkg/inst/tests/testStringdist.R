
library(testthat)

context("General ")
test_that("Argument parsing",{
   expect_equal(stringdist(character(0),"a"),numeric(0))
   expect_equal(stringdist("a",character(0)),numeric(0))
   expect_error(stringdist("a","b",weight=c(-1,1,1,1)))
   expect_error(stringdist("a","b",weight=c(1,0,1,1)))
   expect_error(stringdist("a","b",weight=c(1,1,1,4)))
}) 


context("Optimal String Alignment")
test_that("Edge cases in OSA method",{
   expect_equal(stringdist( "", "",method='osa'),0)
   expect_equal(stringdist( "","a",method='osa'),1)
   expect_equal(stringdist("a", "",method='osa'),1)
   expect_equal(stringdist("a","a",method='osa'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='osa',maxDist=1),-1)
   # Thanks to Daniel Deckhard pointing to this bug
   expect_equal(stringdist("abc","abc",method='osa',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='osa',maxDist=1), -1)
   expect_equal(stringdist("abc","",method='osa',maxDist=2), -1)
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='osa'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='osa'),c(0,1))
})

test_that("weights are handled correctly",{
   # deletion
   expect_equal(stringdist("a","ab", method='osa',weight=c(0.5,1,1,1)),0.5)
   # insertion
   expect_equal(stringdist("ab","a" ,method='osa',weight=c(1,0.5,1,1)),0.5)
   # substitution
   expect_equal(stringdist("b","a" , method='osa',weight=c(1,1,0.5,1)),0.5)
   # transposition
   expect_equal(stringdist("ca","ac",method='osa',weight=c(1,1,1,0.5)),0.5)
   # symmetry property in simple case
   expect_equal(
      stringdist("abc","ac",method='osa',weight=c(0.5,1,1,1)),
      stringdist("ac","abc",method='osa',weight=c(1,0.5,1,1))
   )
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='osa')))
   expect_true(is.na(stringdist('a',NA ,method='osa')))
   expect_true(is.na(stringdist(NA ,NA ,method='osa')))

})

context("Levenstein")
test_that("Edge cases in Levenshtein method",{
   expect_equal(stringdist( "", "",method='lv'),0)
   expect_equal(stringdist( "","a",method='lv'),1)
   expect_equal(stringdist("a", "",method='lv'),1)
   expect_equal(stringdist("a","a",method='lv'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='lv',maxDist=1),-1)
   expect_equal(stringdist("abc","abc",method='lv',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='lv',maxDist=1), -1)
   expect_equal(stringdist("abc","",method='lv',maxDist=2), -1)
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='lv'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='lv'),c(0,1))
})

test_that("weights are handled correctly",{
   # deletion
   expect_equal(stringdist("a","ab", method='lv',weight=c(0.5,1,1)),0.5)
   # insertion
   expect_equal(stringdist("ab","a" ,method='lv',weight=c(1,0.5,1,1)),0.5)
   # substitution
   expect_equal(stringdist("b","a" , method='lv',weight=c(1,1,0.5,1)),0.5)
   # symmetry property in simple case
   expect_equal(
      stringdist("abc","ac",method='lv',weight=c(0.5,1,1,1)),
      stringdist("ac","abc",method='lv',weight=c(1,0.5,1,1))
   )
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='lv')))
   expect_true(is.na(stringdist('a',NA ,method='lv')))
   expect_true(is.na(stringdist(NA ,NA ,method='lv')))
})

context("Damerau-Levenstein")
test_that("Edge cases in DL method",{
   expect_equal(stringdist( "", "",method='dl'),0)
   expect_equal(stringdist( "","a",method='dl'),1)
   expect_equal(stringdist("a", "",method='dl'),1)
   expect_equal(stringdist("a","a",method='dl'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='dl',maxDist=1),-1)
   expect_equal(stringdist("abc","abc",method='dl',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='dl',maxDist=2), -1)
   expect_equal(stringdist("abc","",method='dl',maxDist=2), -1)
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='dl'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='dl'),c(0,1))
})

test_that("weights are handled correctly",{
   # deletion
   expect_equal(stringdist("a","ab", method='dl',weight=c(0.5,1,1,1)),0.5)
   # insertion
   expect_equal(stringdist("ab","a" ,method='dl',weight=c(1,0.5,1,1)),0.5)
   # substitution
   expect_equal(stringdist("b","a" , method='dl',weight=c(1,1,0.5,1)),0.5)
   # transposition
   expect_equal(stringdist("ca","ac",method='dl',weight=c(1,1,1,0.5)),0.5)
   # symmetry property in simple case
   expect_equal(
      stringdist("abc","ac",method='dl',weight=c(0.5,1,1,1)),
      stringdist("ac","abc",method='dl',weight=c(1,0.5,1,1))
   )
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='dl')))
   expect_true(is.na(stringdist('a',NA ,method='dl')))
   expect_true(is.na(stringdist(NA ,NA ,method='dl')))
})

context("Longest Common Substring")
test_that("Edge cases in LCS method",{
   expect_equal(stringdist( "", "",method='lcs'),0)
   expect_equal(stringdist( "","a",method='lcs'),1)
   expect_equal(stringdist("a", "",method='lcs'),1)
   expect_equal(stringdist("a","a",method='lcs'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='lcs',maxDist=1),-1)
   expect_equal(stringdist("abc","abc",method='lcs',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='lcs',maxDist=1), -1)
   expect_equal(stringdist("abc","",method='lcs',maxDist=1), -1)
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='lcs'),c(0,2))
   expect_equal(stringdist('a',c('a','b'),method='lcs'),c(0,2))
})


test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='lcs')))
   expect_true(is.na(stringdist('a',NA ,method='lcs')))
   expect_true(is.na(stringdist(NA ,NA ,method='lcs')))
})


context("Hamming distance")

test_that("Edge cases in DL method",{
   expect_equal(stringdist( "", "",method='h'),0)
   expect_equal(stringdist("a","a",method='h'),0)
})

test_that("Unequal string lengths",{
  expect_equal(stringdist("aa","a",method="h"),-1)
  expect_equal(stringdist("a","aa",method="h"),-1)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='h',maxDist=1),-1)
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='h'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='h'),c(0,1))
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='h')))
   expect_true(is.na(stringdist('a',NA ,method='h')))
   expect_true(is.na(stringdist(NA ,NA ,method='h')))
})

context("Q-gram distance")

test_that("Edge cases in qgram method",{
   expect_equal(stringdist( "", "",method='qgram',q=0), 0)
   expect_equal(stringdist( "", "",method='qgram',q=1),-1)
   expect_equal(stringdist( "","a",method='qgram',q=1),-1)
   expect_equal(stringdist("a", "",method='qgram',q=1),-1)
   expect_equal(stringdist("a","a",method='qgram',q=1), 0)
   expect_error(stringdist("aa","bb",method='qgram',q=-2))
})


test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='qgram',q=1),c(0,2))
   expect_equal(stringdist('a',c('a','b'),method='qgram',q=1),c(0,2))
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='qgram')))
   expect_true(is.na(stringdist('a',NA ,method='qgram')))
   expect_true(is.na(stringdist(NA ,NA ,method='qgram')))
})


test_that("binary tree is cleaned up properly in qgram-tree",{
# explanation: the binary tree storing unique q-grams and q-gram counts is re-used when looping
# over string pairs. (this is not the case with the unsorted lookup table in 'qgram')
  d <- stringdist('abcde',c('edcba','edcba'),method='qgram',q=2)
  expect_equal(d[1],d[2])
})


context("Jaro")
test_that("basic examples and edge cases work",{
  # strings of length 1
  expect_equal(stringdist("a","a",method='jaro'),0);
  expect_equal(stringdist("a","b",method='jaro'),1);
  expect_equal(stringdist("a","",method='jaro'), 1);
  expect_equal(stringdist("","",method='jaro'), 0);
})

test_that("Extended examples work",{
  # cases from wikipedia
  expect_equal(
    round(1-stringdist("martha","marhta",method='jaro'),3),
    0.944
  )
  expect_equal(
    round(1-stringdist("dixon","dicksonx",method='jaro'),3),
    0.767
  )
  expect_equal(
    round(1-stringdist("duane","dwayne",method='jaro'),3),
    0.822
  )
  # Other cases
  # 4 matches, no transpositions, short first string with non-matching character.
  expect_equal(stringdist("axiou","aaeeiioouu",method='jaro'),1-(4/5+4/10 + 4/4)/3);
  # non-matching characters in both strings
  expect_equal(stringdist("abcdeu","abxde",method='jaro'),1-(4/6+4/5+4/4)/3);
})

test_that("distance is symmetric",{
  expect_equal(
    round(stringdist("martha","marhta",method='jaro'),8),
    round(stringdist("marhta","martha",method='jaro'),8)
  )
  expect_equal(
    round(stringdist("dicksonx","dixon",method='jaro'),8),
    round(stringdist("dixon","dicksonx",method='jaro'),8)
  )
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='jaro'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='jaro'),c(0,1))
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='jaro')))
   expect_true(is.na(stringdist('a',NA ,method='jaro')))
   expect_true(is.na(stringdist(NA ,NA ,method='jaro')))
})


context("stringdistmatrix")
test_that("dimensions work out",{
    expect_equivalent(
        dim(stringdistmatrix(c("aa","bb","cc"),c("aa","cc"))),
        c(3,2)
    )
})


