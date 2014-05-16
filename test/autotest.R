library(testthat)


dyn.load("../pkg/src/stringdist.so")
auto_test("../pkg/R", "../pkg/inst/tests")

