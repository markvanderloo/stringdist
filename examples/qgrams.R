
qgrams('hello world',q=3)

# q-grams are counted uniquely over a character vector
qgrams(rep('hello world',2),q=3)

# to count them separately, do something like
x <- c('hello', 'world')
lapply(x,qgrams, q=3)

<<<<<<< HEAD
# a tonque twister
x <- "peter piper picked a peck of pickled peppers"
qgrams(x, q=2, useNames=FALSE) 
qgrams(x, q=2, useNames=FALSE) 
qgrams(x, q=2, useBytes=TRUE)
qgrams(x, q=2, useBytes=TRUE, useNames=TRUE)
=======
# output rows may be named, and you can pass any number of character vectors
x <- "I will not buy this record, it is scratched"
y <- "My hovercraft is full of eels"
z <- c("this", "is", "a", "dead","parrot")
qgrams(A = x, B = y, C = z,q=2)

>>>>>>> 9483d4021b21dee5286f2b8127186a4d7df97209


