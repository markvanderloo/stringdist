
library(testthat)

### -------------------------------------------------------------
context("General ")
test_that("Argument parsing",{
   expect_equal(stringdist(character(0),"a"),numeric(0))
   expect_equal(stringdist("a",character(0)),numeric(0))
   expect_error(stringdist("a","b",weight=c(-1,1,1,1)))
   expect_error(stringdist("a","b",weight=c(1,0,1,1)))
   expect_error(stringdist("a","b",weight=c(1,1,1,4)))
   expect_warning(stringdist(letters[1:3],letters[1:2]))
   expect_warning(stringdist(list('a'),'a'))
   expect_warning(stringdist('a',list('a')))
   expect_warning(stringdistmatrix(list('a')))
   expect_warning(stringdistmatrix(list('a'),list('b')))
}) 


### -------------------------------------------------------------
context("Optimal String Alignment")
test_that("Edge cases in OSA method",{
   expect_equal(stringdist( "", "",method='osa'),0)
   expect_equal(stringdist( "","a",method='osa'),1)
   expect_equal(stringdist("a", "",method='osa'),1)
   expect_equal(stringdist("a","a",method='osa'),0)
   # Thanks to Frank Binder for poining out this bug
   expect_equal(stringdist("ab","aba",method='osa'),1)
   # Thanks to Lauri Koobas for pointing out this bug.
   expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd")))),1)


})

test_that("max distance yields warning",{
  expect_warning(stringdist("abc","abc",method='osa',maxDist=1))
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
   
  expect_equal( # two deletions from source (b) to target (a)
    stringdist("","aa",weight=c(0.5,1,1,1),method="osa"),1
  )
  expect_equal( # two deletions from source (b) to target (a)
    stringdist("","aa",weight=c(0.5,1,1,1),method="lv"),1
  )
  expect_equal( # two deletions from source (b) to target (a)
    stringdist("","aa",weight=c(0.5,1,1,1),method="dl"),1
  )
  
  # Thanks to Zach Price for reporting this bug.
  expect_equal(
    stringdist("ABC", "BC", method = "lv", weight = c(i=.1, d=.1, s=.1)),.1
  )
  expect_equal(
    stringdist("ABC", "BC", method = "lv", weight = c(i=.1, d=.1, s=1)),.1
  )
  
  expect_equal(
    stringdist("ABC", "BC", method = "osa", weight = c(i=.1, d=.1, s=.1,t=.1)),.1
  )
  expect_equal(
    stringdist("ABC", "BC", method = "osa", weight = c(i=.1, d=.1, s=1,t=.1)),.1
  )
  expect_equal(
    stringdist("ABC", "BC", method = "dl", weight = c(i=.1, d=.1, s=.1,t=.1)),.1
  )
  expect_equal(
    stringdist("ABC", "BC", method = "dl", weight = c(i=.1, d=.1, s=1,t=.1)),.1
  )
  # examples from the paper; Tanks to Nathalia Potocka for reporting.
  expect_equal(stringdist("leia","leela",method="lv",weight=c(i=.1,d=1,s=1)),1.1)
  expect_equal(stringdist("leia","leela",method="lv",weight=c(i=1,d=.1,s=1)),2)
  expect_equal(stringdist("a","b",method="lv",weight=c(i=.1,d=1,s=.3)),.3)
  expect_equal(stringdist("a","b",method="osa",weight=c(i=.1,d=1,s=.3,1)),.3)
  expect_equal(stringdist("a","b",method="dl",weight=c(i=.1,d=1,s=.3,t=1)),.3)
  expect_equal(stringdist("leia","leela",method="dl",weight=c(i=1,d=.1,s=1,t=1)),2)

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
   expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="lv"))),1)
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
   expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="lcs"))),1)
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
   expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="h"))),1)
})

test_that("Unequal string lengths",{
  expect_equal(stringdist("aa","a",method="h"),Inf)
  expect_equal(stringdist("a","aa",method="h"),Inf)
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
   expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="qgram"))),1)
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
  # numerical accuracy test (thanks to Ben Haller)
  # note that 1 - 2/(sqrt(2)*sqrt(2)) != 0, so this used to give ~2.2E-16. 
  expect_equal( stringdist("ab","ab",method="cosine"),0.0,tolerance=0.0 ) 
  expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="cosine"))),1)
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
  expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="jaccard"))),1)
})


### -------------------------------------------------------------
context("Jaro")
test_that("basic examples and edge cases work",{
  # strings of length 1
  expect_equal(stringdist("a","a",method='jw'),0);
  expect_equal(stringdist("a","b",method='jw'),1);
  expect_equal(stringdist("a","",method='jw'), 1);
  expect_equal(stringdist("","",method='jw'), 0);
  # following test added after a bug report of Carol Gan:
  expect_equal(stringdist("tire","tree",method="jw"),stringdist("tree","tire",method="jw"));
  expect_equal(sum(is.na(stringdist(c("a", NA, "b", "c"), c("aa", "bb", "cc", "dd"),method="jw"))),1)
  # following test added after issue #42 by github user gtumuluri
  # thanks to Jan for providing this simple example.
  expect_equal(stringdist("DHCXXXXX","HCDXXXXX",method="jw"),stringdist("HCDXXXXX","DHCXXXXX",method="jw"))
  # the following test added after issue #42, comment by desource90
  expect_equal(
    stringdist("RICK WARREN","WARREN BUFFET",method="jw")
    , 1 - (1/3)*(7/13 + 7/11 + (7-3.5)/7))
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
  # idem, with weights
  expect_equal(
    round(1 - stringdist("crate","trace",metho='jw',weight=c(0.5,1,1)),8),
    round((0.5*3/5 + 3/5 + (3-0)/3)/3,8)
  )
  expect_equal(
    round(1 - stringdist("crate","trace",metho='jw',weight=c(1,0.5,1)),8),
    round((3/5 + 0.5*3/5 + (3-0)/3)/3,8)
  )
  expect_equal(
    round(1 - stringdist("crate","trace",metho='jw',weight=c(1,1,0.5)),8),
    round((3/5 + 3/5 + 0.5*(3-0)/3)/3,8)
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
    expect_equivalent(
        dim(stringdistmatrix(c("aa","bb","cc"),c("aa","cc"),useBytes=TRUE)),
        c(3,2)
    ) 
    expect_equivalent( # bug #28
      dim(stringdistmatrix('foo',letters[1:3])), c(1,3)
    )
})

test_that("stringdistmatrix-lower-tri can output long vectors",{
   # skipped on CRAN because of high memory use.
   skip_on_cran()
   # Error when input vector yields a vector too big for a long vector.
   out <- tryCatch(stringdistmatrix(character(100663296+1),method="hamming")
      , error = function(e) e$message )
   expect_equal(class(out),"character")
   expect_match(out, "exceeds maximum allowed")
})

test_that('stringdistmatrix yields correct distances',{
  x <- paste0('Mot',intToUtf8(0x00F6),'rhead') # correct spelling
  y <- 'Motorhead' # Pissing off Lemmy.
  v <- c(x,y)
  d11 <- stringdist(x,x)
  d12 <- stringdist(x,y)
  d22 <- stringdist(y,y)
  expect_equal(
    stringdistmatrix(v,v)
    , matrix(c(d11,d12,d12,d22),nrow=2,ncol=2)
  )
  d11 <- stringdist(x,x,useBytes=TRUE)
  d12 <- stringdist(x,y,useBytes=TRUE)
  d22 <- stringdist(y,y,useBytes=TRUE)
  expect_equal(
    stringdistmatrix(v,v,useBytes=TRUE)
    , matrix(c(d11,d12,d12,d22),nrow=2,ncol=2)
  )
})

test_that("stringdistmatrix gives correct labels",{
  a <- c(k1="jan",k2="pier",k3="joris")
  b <- c(f1 = "jip", f2="janneke")
  expect_equal(
    dimnames(stringdistmatrix(a,b,useNames=TRUE))
    , list(as.character(a),as.character(b))
  )
  expect_equal(
    dimnames(stringdistmatrix(a,b,useNames="strings"))
    , list(as.character(a),as.character(b))
  )  
  expect_equal(
    dimnames(stringdistmatrix(a,b,useNames="names"))
    , list(c("k1","k2","k3"),c("f1","f2"))
  )
  
})


test_that("stringdistmatrix with single argument",{
  d <- stringdistmatrix(c("aap","noot","mies","boom","roos","vis"))
  expect_equal(class(d),"dist")
  expect_equal(length(d),15)
  d <- stringdistmatrix(c("a",NA,"b"))
  expect_equal(sum(is.na(d)),2)
  
  a <- c(k1 = "aap",k2="noot")
  expect_true(is.null(attr(stringdistmatrix(a,useNames="none"),"Labels")))
  expect_identical(stringdistmatrix(a,useNames="none"),stringdistmatrix(a,useNames=FALSE))
  expect_equal(
      attr(stringdistmatrix(a,useNames="strings"),"Labels")
    , as.character(a)
  )

  expect_equal(
    attr(stringdistmatrix(a,useNames=TRUE), "Labels")
    , c("aap","noot")
  )
  
  expect_equal(
    attr(stringdistmatrix(a,useNames="names"), "Labels")
    , c("k1","k2")
  )
  
  
})


context("stringdist: useBytes")
test_that("useBytes gets NA",{
  expect_true(is.na(stringdist('a',NA,method='osa',useBytes=TRUE)))
  expect_true(is.na(stringdist('a',NA,method='lv',useBytes=TRUE)))
  expect_true(is.na(stringdist('a',NA,method='dl',useBytes=TRUE)))
  expect_true(is.na(stringdist('a',NA,method='hamming',useBytes=TRUE)))
})

test_that("useBytes translates correctly to numeric",{
  # smoketest
  set.seed(1)
  x <- sapply(sample(5:25,10,replace=TRUE),function(x) paste(letters[x],collapse=""))
  y <- sample(x)
  expect_equal(
    stringdist(x,y,method='osa',useBytes=TRUE)
  , stringdist(x,y,method='osa',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='lv',useBytes=TRUE)
  , stringdist(x,y,method='lv',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='dl',useBytes=TRUE)
  , stringdist(x,y,method='dl',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='hamming',useBytes=TRUE)
  , stringdist(x,y,method='hamming',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='jw',useBytes=TRUE)
  , stringdist(x,y,method='jw',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='lcs',useBytes=TRUE)
  , stringdist(x,y,method='lcs',useBytes=FALSE))
  expect_equal(
    stringdist(x,y,method='qgram',q=3,useBytes=TRUE)
  , stringdist(x,y,method='qgram',q=3,useBytes=FALSE))

})

test_that("useBytes really analyses bytes",{
  x <- paste0('Mot',intToUtf8(0x00F6),'rhead') # correct spelling
  y <- 'Motorhead' # Pissing off Lemmy.
  expect_equal(stringdist(x,y,method='dl',useBytes=TRUE),  2)
  expect_equal(stringdist(x,y,method='hamming',useBytes=TRUE),  Inf)
  expect_equal(stringdist(x,y,method='osa',useBytes=TRUE), 2)
  expect_equal(stringdist(x,y,method='lv',useBytes=TRUE),  2)
  expect_equal(
    round(stringdist(x,y,method='jw',useBytes=TRUE),3),
    round(1-(1/3)*(8/9 + 8/10 + 1),3)
  )
  expect_equal(stringdist(x,y,method='lcs',useBytes=TRUE),  3)
  expect_equal(stringdist(x,y,method='qgram',q=3,useBytes=TRUE),  7)
})



### -------------------------------------------------------------
context("Soundex distance")

test_that("",{
  expect_equal(stringdist("", "0000",method='soundex'),0)
  expect_equal(stringdist("john","jan",method='soundex'),0)
  expect_equal(stringdist("schoen","son",method='soundex'),0)
  expect_equal(stringdist("ssssss","sa",method='soundex'),0)
  expect_equal(stringdist("bfpv","ba",method='soundex'),0)
  # cgjkqsxz receive same code
  expect_equal(stringdist("cgjkqsxz","ca",method='soundex'),0)
  # dt receive same code
  expect_equal(stringdist("dt","da",method='soundex'),0)  
  # mn receive same code; l has seperate code
  expect_equal(stringdist("lmn","lam",method='soundex'),0)  
  # h and w are ignored
  expect_equal(stringdist("rhw","r",method='soundex'),0)  
  # vowels are removed
  expect_equal(stringdist("raeiouym","rm",method='soundex'),0)  
  # non-letters are ignores
  expect_equal(stringdist("r00d","rt",method='soundex'),0)  
  # consonants are not merged when a vowel in between
  expect_equal(stringdist("sock", "sck", method='soundex'),1)
  x <- "Motorhead"
  y <- paste0("Mot",intToUtf8(0x00F6),"rhead") # with o-umlaut
  expect_warning(stringdist(x,y,method='soundex',useBytes=TRUE))  
})

test_that("Shortest argument is recycled",{
  expect_equal(stringdist(c('a','b'),'a',method='soundex'),c(0,1))
  expect_equal(stringdist('a',c('a','b'),method='soundex'),c(0,1))
})

test_that("NA's are handled correctly",{
  expect_true(is.na(stringdist(NA ,'a',method='soundex')))
  expect_true(is.na(stringdist('a',NA ,method='soundex')))
  expect_true(is.na(stringdist(NA ,NA ,method='soundex')))
})

test_that("non-printable ascii and non-ascii encoding is detected",{
  ouml <- intToUtf8("0x00F6")
  # non-ascii within string
  x <- paste0("Mot",ouml,"rhead")
  y <- paste0(ouml,"zzy")
  # non-printable ascii's in string
  z <- paste0("\r","hello","\t")
  # business as usual
  expect_equal(stringdist('Ozzy','Lemmy',method='soundex'),1)
  # first argument triggers warning
  expect_warning(stringdist(x,'Lemmy',method='soundex'))
  expect_warning(stringdist(y,'Ozzy',method='soundex'))
  expect_warning(stringdist(z,'Ozzy',method='soundex'))
  
  # second argument triggers warning
  expect_warning(stringdist('Lemmy',x,method='soundex'))
  expect_warning(stringdist('Ozzy',y,method='soundex'))
  expect_warning(stringdist('Ozzy',z,method='soundex'))
})



