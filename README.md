stringdist
==========

String distance functions for R

String distance functions are scattered around R, and R's packages. This package
is a re-implementation of some common (weighted) distance functions, offered
through a uniform interface. So far distance functions include:

* Hamming distance; 
* Levenshtein distance (weighted);
* Restricted Damerau-Levenshtein distance (weighted);
* Full Damerau-Levenshtein distance (weighted).

To my best knowledge, the latter was not available before in R. Workhorse functions are implemented in C. The package offers two main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix between two input character vectors, optionally running in parallel.

TODO
----
* distances based on ngrams
* jaro-winkler distance
* ~~optionally use user-defined cluster for parallel computations~~
* small performance tweeks, like detect where |nchar(a)-nchar(b)| > maxDist or smarter distribution of jobs over clusters

