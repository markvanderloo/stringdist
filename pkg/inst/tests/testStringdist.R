
library(testthat)

### -------------------------------------------------------------
context("General ")
test_that("Argument parsing",{
   expect_equal(stringdist(character(0),"a"),numeric(0))
   expect_equal(stringdist("a",character(0)),numeric(0))
   expect_error(stringdist("a","b",weight=c(-1,1,1,1)))
   expect_error(stringdist("a","b",weight=c(1,0,1,1)))
   expect_error(stringdist("a","b",weight=c(1,1,1,4)))
}) 


### -------------------------------------------------------------
context("Optimal String Alignment")
test_that("Edge cases in OSA method",{
   expect_equal(stringdist( "", "",method='osa'),0)
   expect_equal(stringdist( "","a",method='osa'),1)
   expect_equal(stringdist("a", "",method='osa'),1)
   expect_equal(stringdist("a","a",method='osa'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='osa',maxDist=1), Inf)
   # Thanks to Daniel Deckhard pointing to this bug
   expect_equal(stringdist("abc","abc",method='osa',maxDist=1), 0)
   expect_equal(stringdist("aa","bbb",method='osa',maxDist=2), Inf)
   expect_equal(stringdist("bbb","aa",method='osa',maxDist=2), Inf)
   expect_equal(stringdist("","abc",method='osa',maxDist=1), Inf)
   expect_equal(stringdist("abc","",method='osa',maxDist=2), Inf)
})

test_that("transpositions are found",{
  expect_equal(stringdist("ab","ba",method='osa'),1)
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

### -------------------------------------------------------------
context("Levenstein")
test_that("Edge cases in Levenshtein method",{
   expect_equal(stringdist( "", "",method='lv'),0)
   expect_equal(stringdist( "","a",method='lv'),1)
   expect_equal(stringdist("a", "",method='lv'),1)
   expect_equal(stringdist("a","a",method='lv'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='lv',maxDist=1), Inf)
   expect_equal(stringdist("abc","abc",method='lv',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='lv',maxDist=1), Inf)
   expect_equal(stringdist("abc","",method='lv',maxDist=2), Inf)
   expect_equal(stringdist("aa","bbb",method='lv',maxDist=2), Inf)
   expect_equal(stringdist("bbb","aa",method='lv',maxDist=2), Inf)
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

### -------------------------------------------------------------
context("Damerau-Levenstein")
test_that("Edge cases in DL method",{
   expect_equal(stringdist( "", "",method='dl'),0)
   expect_equal(stringdist( "","a",method='dl'),1)
   expect_equal(stringdist("a", "",method='dl'),1)
   expect_equal(stringdist("a","a",method='dl'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='dl',maxDist=1),Inf)
   expect_equal(stringdist("abc","abc",method='dl',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='dl',maxDist=2), Inf)
   expect_equal(stringdist("abc","",method='dl',maxDist=2), Inf)
   expect_equal(stringdist("aa","bbb",method='dl',maxDist=2), Inf)
   expect_equal(stringdist("bbb","aa",method='dl',maxDist=2), Inf)
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

### -------------------------------------------------------------
context("Longest Common Substring")
test_that("Edge cases in LCS method",{
   expect_equal(stringdist( "", "",method='lcs'),0)
   expect_equal(stringdist( "","a",method='lcs'),1)
   expect_equal(stringdist("a", "",method='lcs'),1)
   expect_equal(stringdist("a","a",method='lcs'),0)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='lcs',maxDist=1),Inf)
   expect_equal(stringdist("abc","abc",method='lcs',maxDist=1), 0)
   expect_equal(stringdist("","abc",method='lcs',maxDist=1), Inf)
   expect_equal(stringdist("abc","",method='lcs',maxDist=1), Inf)
   expect_equal(stringdist("aa","bbb",method='lcs',maxDist=2), Inf)
   expect_equal(stringdist("bbb","aa",method='lcs',maxDist=2), Inf)
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


### -------------------------------------------------------------
context("Hamming distance")
test_that("Edge cases in DL method",{
   expect_equal(stringdist( "", "",method='h'),0)
   expect_equal(stringdist("a","a",method='h'),0)
})

test_that("Unequal string lengths",{
  expect_equal(stringdist("aa","a",method="h"),Inf)
  expect_equal(stringdist("a","aa",method="h"),Inf)
})

test_that("max distance is obeyed",{
   expect_equal(stringdist("aa","bb",method='h',maxDist=1),Inf)
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



### -------------------------------------------------------------
context("Q-gram distance")

test_that("Edge cases in qgram method",{
   expect_equal(stringdist( "", "",method='qgram',q=0), 0)
   expect_equal(stringdist( "", "",method='qgram',q=1),Inf)
   expect_equal(stringdist( "","a",method='qgram',q=1),Inf)
   expect_equal(stringdist("a", "",method='qgram',q=1),Inf)
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

### -------------------------------------------------------------
context("cosine distance") 
# basic engine is q-gram so we need limited testing
test_that("cosine distance computes correctly",{
  expect_equal(
    round(stringdist("aaa","abc",method="cosine",q=1),8),
    round(1-1/sqrt(3),8)
  )
  expect_equal(
    round(stringdist("aaa","abc",method="cosine",q=2),8),
    1.0
  )
})

context("Jaccard distance")
# basic engine is q-gram so we need limited testing
test_that("Jaccard distance computes correctly",{
  expect_equal(
    round(stringdist("aaa","abc",method="jaccard",q=1),8),
    round(1-1/3,8)
  )
  expect_equal(
    round(stringdist("aaa","abc",method="jaccard",q=2),8),
    1.0
  )
})


### -------------------------------------------------------------
context("Jaro")
test_that("basic examples and edge cases work",{
  # strings of length 1
  expect_equal(stringdist("a","a",method='jw'),0);
  expect_equal(stringdist("a","b",method='jw'),1);
  expect_equal(stringdist("a","",method='jw'), 1);
  expect_equal(stringdist("","",method='jw'), 0);
})

test_that("Extended examples work",{
  # cases from wikipedia
  expect_equal(
    round(1-stringdist("martha","marhta",method='jw'),3),
    0.944
  )
  expect_equal(
    round(1-stringdist("dixon","dicksonx",method='jw'),3),
    0.767
  )
  expect_equal(
    round(1-stringdist("duane","dwayne",method='jw'),3),
    0.822
  )
  # out-of range matches
  expect_equal(
    round(1 - stringdist("crate","trace",metho='jw'),8),
    round((3/5 + 3/5 + (3-0)/3)/3,8)
  )

  # Other cases
  # 4 matches, no transpositions, short first string with non-matching character.
  expect_equal(stringdist("axiou","aaeeiioouu",method='jw'),1-(4/5+4/10 + 4/4)/3);
  # non-matching characters in both strings
  expect_equal(stringdist("abcdeu","abxde",method='jw'),1-(4/6+4/5+4/4)/3);
})

test_that("distance is symmetric",{
  expect_equal(
    round(stringdist("martha","marhta",method='jw'),8),
    round(stringdist("marhta","martha",method='jw'),8)
  )
  expect_equal(
    round(stringdist("dicksonx","dixon",method='jw'),8),
    round(stringdist("dixon","dicksonx",method='jw'),8)
  )
})

test_that("Shortest argument is recycled",{
   expect_equal(stringdist(c('a','b'),'a',method='jw'),c(0,1))
   expect_equal(stringdist('a',c('a','b'),method='jw'),c(0,1))
})

test_that("NA's are handled correctly",{
   expect_true(is.na(stringdist(NA ,'a',method='jw')))
   expect_true(is.na(stringdist('a',NA ,method='jw')))
   expect_true(is.na(stringdist(NA ,NA ,method='jw')))
})

### -------------------------------------------------------------
context("Jaro-Winkler")
test_that("wikipedia examples",{
  expect_equal(
    round(stringdist("martha","marhta",method="jw",p=0.1),3),
    1-0.961
  )
  expect_equal(
    round(stringdist("dwayne","duane",method="jw",p=0.1),3),
    round(1-0.84,3)
  )
  expect_equal(
    round(stringdist("dixon","dicksonx",method="jw",p=0.1),3),
    round(1-0.813,3)
  )
})


context("stringdistmatrix")
test_that("dimensions work out",{
    expect_equivalent(
        dim(stringdistmatrix(c("aa","bb","cc"),c("aa","cc"))),
        c(3,2)
    )
})

# since the result of useBytes depends on the encoding used we cannot really unit-test on
# distance values. However, we can do some basic tests and check for crashes.
context("useBytes")
test_that("useBytes gets NA",{
  expect_true(is.na(stringdist('a',NA,method='osa',useBytes=TRUE)))
  expect_true(is.na(stringdist('a',NA,method='lv',useBytes=TRUE)))
})

test_that("useBytes doesn't crash",{
  # smoketest
  x <- sapply(sample(5:25,10,replace=TRUE),function(x) paste(letters[x],collapse=""))
  stringdist(x,sample(x),useBytes=TRUE)

})
