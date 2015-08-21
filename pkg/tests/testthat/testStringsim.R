library(testthat)

context("stringsim")

# We expect that two completely different strings have a similarity of
# 0 and two completely equal strings a similarity of 1
methods <- c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")
for (method in methods) {

  test_that(paste0("Similarity is between 0 and 1 for ", method), {
    expect_that(stringsim("bb", "cc", method=method), equals(0))
    expect_that(stringsim("bb", "bb", method=method), equals(1))
  })
}

for (method in methods[c(1:5,9:10)]){
  test_that(paste0("Edge cases for ", method), {
    expect_that(stringsim(c("a", ""), "", method=method), equals(c(0, 1)))
    
    expect_that(stringsim(c("kkk", "bbb"), "bbb", method=method), 
      equals(stringsim("bbb", c("kkk", "bbb"), method=method)))
  })
}
for (method in methods[6:8]){
  test_that(paste0("Edge cases for ", method), {
    expect_that(stringsim(c("a", ""), "", method=method,q=0), equals(c(1, 1)))
    
    expect_that(stringsim(c("kkk", "bbb"), "bbb", method=method,q=0), 
                equals(stringsim("bbb", c("kkk", "bbb"), method=method,q=0)))
  })
}




test_that("Stringsim gets correct values with or without useBytes",{
  x <- "ao"
  y <- paste0("a",intToUtf8(0x00F6)) # o-umlaut
  expect_equal(stringsim(x,y,method="osa", useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="osa", useBytes=TRUE ), 1-2/3)
  expect_equal(stringsim(x,y,method="lv", useBytes=FALSE),  1-1/2)
  expect_equal(stringsim(x,y,method="lv", useBytes=TRUE ),  1-2/3)
  expect_equal(stringsim(x,y,method="dl", useBytes=FALSE),  1-1/2)
  expect_equal(stringsim(x,y,method="dl", useBytes=TRUE ),  1-2/3)
  expect_equal(stringsim(x,y,method="hamming", useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="hamming", useBytes=TRUE ), 1-1)
  expect_equal(stringsim(x,y,method="lcs", useBytes=FALSE),     1-1/2)
  expect_equal(stringsim(x,y,method="lcs", useBytes=TRUE ),     1-3/5)
  expect_equal(stringsim(x,y,method="qgram", q=1, useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="qgram", q=1, useBytes=TRUE ), 1-3/5)
  expect_equal(stringsim(x,y,method="cosine", q=1, useBytes=FALSE), 1-1/2)
  expect_equal(stringsim(x,y,method="cosine", q=1, useBytes=TRUE ), 1-(1-1/sqrt(6)))
  expect_equal(stringsim(x,y,method="jaccard", q=1, useBytes=FALSE), 1-2/3)
  expect_equal(stringsim(x,y,method="jaccard", q=1, useBytes=TRUE ), 1-3/4)
  expect_equal(stringsim(x,y,method="jw", useBytes=FALSE), 1-1/3)
  expect_equal(stringsim(x,y,method="jw", useBytes=TRUE ), (1/2 + 1/3 +1)/3)
})  

context("seq_sim")

test_that("elementary seq_sim test",{
  expect_equal(
    seq_sim(list(1:3,2:4),list(1:3))
     , stringsim(c("abc","bcd"),"abc") 
  )
})




