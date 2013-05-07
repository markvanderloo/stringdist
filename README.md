stringdist
==========

String distance functions for R

String distance functions are scattered around R, and R's packages. This package
is a re-implementation of some common (weighted) distance functions, offered
through a uniform interface. As of version 0.5.0, distance functions include:

* Hamming distance; 
* Levenshtein distance (weighted);
* Restricted Damerau-Levenshtein distance (weighted, a.k.a. Optimal String Alignment);
* Full Damerau-Levenshtein distance (weighted);
* Longest Common Substring distance;
* Q-gram distance (two implementations).

To my best knowledge, the latter three were not available before in R. Workhorse
functions are implemented in C. The package offers two main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix between two input character vectors, optionally running in parallel.

TODO
----
* Episode distance
* ~~Longest common substring~~
* ~~distances based on q-grams~~ 
* jaro-winkler distance
* ~~optionally use user-defined cluster for parallel computations~~
* small performance tweeks, like detect where |nchar(a)-nchar(b)| > maxDist or smarter distribution of jobs over clusters

Could
----
Perhaps in the future I'll add supporting functionality such as
- An encoding sniffer, detecting character encodings from files (not sure if that's avaiable in R already)
- Some string normalizing functionality
- ...

