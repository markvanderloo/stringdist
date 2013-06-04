
get_qgrams('hello world',q=3)

# q-grams are counted uniquely over a character vector
get_qgrams(rep('hello world',2),q=3)

# to count them separately, do something like
x <- c('hello', 'world')
lapply(x,get_qgrams, q=3)




