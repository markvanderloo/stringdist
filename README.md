
[![Build Status](https://travis-ci.org/markvanderloo/stringdist.svg?branch=master)](https://travis-ci.org/markvanderloo/stringdist)
[![Coverage Status](https://coveralls.io/repos/markvanderloo/stringdist/badge.svg)](https://coveralls.io/r/markvanderloo/stringdist) 
[![CRAN](http://www.r-pkg.org/badges/version/stringdist)](http://cran.r-project.org/web/packages/stringdist/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/stringdist)](http://cran.r-project.org/package=stringdist/) 


## stringdist

* Approximate matching and string distance calculations for R. 
* All distance and matching operations are system- and encoding-independent.

The package offers the following main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix for one or two vectors
* `stringsim` computes a string similarity between 0 and 1, based on `stringdist`
* `amatch` is a fuzzy matching equivalent of R's native `match` function
* `ain` is a fuzzy matching equivalent of R's native `%in%` operator

These functions are built upon `C`-code that re-implements some common (weighted) string
distance functions. Distance functions include:

* Hamming distance; 
* Levenshtein distance (weighted);
* Restricted Damerau-Levenshtein distance (weighted, a.k.a. Optimal String Alignment);
* Full Damerau-Levenshtein distance (weighted);
* Longest Common Substring distance;
* Q-gram distance
* cosine distance for q-gram count vectors (= 1-cosine similarity)
* Jaccard distance for q-gram count vectors (= 1-Jaccard similarity)
* Jaro, and Jaro-Winker distance
* Soundex-based string distance.

Also, there are some utility functions:

* `qgrams()` tabulates the qgrams in one or more `character` vectors.
* `phonetic()` computes phonetic codes of strings (currently only soundex)
* `printable_ascii()` is a utility function that detects non-printable ascii or non-ascii characters.



#### Installation

To install the latest release from CRAN, open an R terminal and type

`install.packages('stringdist')`

Beta versions are released through my [drat](http://www.r-pkg.org/pkg/drat) repository. These versions build and pass all current tests correctly on Linux
but builds have not been tested on all architectures that [CRAN](http://cran.r-projecet.org) supports. Windows users will also need to
have [rtools](http://cran.r-project.org/bin/windows/Rtools/) installed.

```
drat::addRepo("markvanderloo")
install.packages("stringdist")
```

To obtain the package from the very latest source code open a `bash` terminal (or `git bash` if you work under Windows
with `msysgit`) and type

```
git clone https://github.com/markvanderloo/stringdist.git
cd stringdist
bash ./build.bash
R CMD INSTALL output/stringdist_*.tar.gz
```

Warning: the github version can change any time and may not even build properly. As most
of the code is written in `C`, the development version may crash your `R`-session.



#### Resources

* A [paper](http://journal.r-project.org/archive/2014-1/loo.pdf) on stringdist has been published in the R-journal * [Slides](http://www.slideshare.net/MarkVanDerLoo/stringdist-use-r2014) of te _useR!2014_ conference.

#### Note to users: deprecated arguments as of >= 0.9.0, >= 0.9.2

Parallelization used to be based on R's ```parallel``` package, that works by spawning several R sessions in the background. As of version 0.9.0, ```stringdist``` uses the more efficient ```openMP``` protocol to parallelize everything under the hood. 

The following arguments have become obsolete and will be removed somewhere in 2016:
* Argument `cluster` for function `stringdistmatrix`.
* Argument `maxDist` for functions `stringdist` and `stringdistmatrix` (not `amatch`).
* Argument `ncores` for function `stringdistmatrix` 


