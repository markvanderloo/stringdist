options(sd_num_thread=1L)

x <- c("ca", "abc", "cba")
expect_equal(stringsimmatrix(x), t(stringsimmatrix(x)))



 
