# this would crash R because of over-asking memory
# it depends on the system really, so we only run this at the
# comfort of our home
if (FALSE){
  x <- paste(letters[sample(1:length(letters),32800,replace=TRUE)], collapse="")
  expect_error(stringdist(x,x))
}


