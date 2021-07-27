
# x <- "IÑIGO", we avoid problems on Windows here.
x <- intToUtf8(c(73, 209,  73,  71,  79))

expect_equal(stringdist("INIGO", x, method="lv", useBytes=FALSE),1)
expect_equal(amatch("INIGO", x, method="lv",maxDist=1),1)



