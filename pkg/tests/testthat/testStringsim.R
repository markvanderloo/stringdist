library(testthat)

context("stringsim")

# We expect that two completely different strings have a similarity of
# 0 and two completely equal strings a similarity of 1
for (method in c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", 
    "jaccard", "jw", "soundex")) {

  test_that(paste0("Similarity is between 0 and 1 for ", method), {
    expect_that(stringsim("bb", "cc", method=method), equals(0))
    expect_that(stringsim("bb", "bb", method=method), equals(1))
  })

  test_that(paste0("Edge cases for ", method), {
    expect_that(stringsim(c("a", ""), "", method=method, q=0), equals(c(0, 1)))
    
    expect_that(stringsim(c("kkk", "bbb"), "bbb", method=method, q=0), 
      equals(stringsim("bbb", c("kkk", "bbb"), method=method, q=0)))
  })
}




