stringdist
==========

String distance functions for R

String distance functions are scattered around R, and R's packages. This package
is a re-implementation of some common (weighted) distance functions, offered
through a uniform interface. As of version 0.5.0, distance functions include:

* Hamming distance; 
* Levenshtein distance (weighted);
* Restricted Damerau-Levenshtein distance (weighted);
* Full Damerau-Levenshtein distance (weighted);
* Q-gram distance.

To my best knowledge, the latter two were not available before in R. Workhorse
functions are implemented in C. The package offers two main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix between two input character vectors, optionally running in parallel.

TODO
----
* Episode distance
* Longest common subsequence
* ~~distances based on q-grams~~ A simple implementation of the Q-gram distance is now present.
* jaro-winkler distance
* ~~optionally use user-defined cluster for parallel computations~~
* small performance tweeks, like detect where |nchar(a)-nchar(b)| > maxDist or smarter distribution of jobs over clusters

