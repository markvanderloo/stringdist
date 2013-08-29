
qgrams('hello world',q=3)

# q-grams are counted uniquely over a character vector
qgrams(rep('hello world',2),q=3)

# to count them separately, do something like
x <- c('hello', 'world')
lapply(x,qgrams, q=3)

# a tonque twister
x <- "peter piper picked a peck of pickled peppers"
qgrams(x, q=2, useNames=FALSE) 
qgrams(x, q=2, useNames=FALSE) 
qgrams(x, q=2, useBytes=TRUE)
qgrams(x, q=2, useBytes=TRUE, useNames=TRUE)


