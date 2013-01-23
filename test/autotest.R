library(testthat)



dyn.load("../pkg/src/dl.so")
auto_test("../pkg/R", "../pkg/inst/tests")

