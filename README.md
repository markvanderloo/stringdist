stringdist
==========

* Approximate matching and string distance calculations for R. 
* All distance and matching operations are system- and encoding-independent.

Current CRAN version: **0.8.0** [NEWS](http://cran.r-project.org/web/packages/stringdist/NEWS)

The package offers four main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix between two input character vectors, optionally running in parallel.
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

Besides the above the functions, `qgrams` tabulates the qgrams in a `character` vector.

Resources
----------
* A [paper](http://journal.r-project.org/archive/2014-1/loo.pdf) on stringdist has been published in the R-journal 
* [Slides](http://www.slideshare.net/MarkVanDerLoo/stringdist-use-r2014) of te _useR!2014_ conference.


Note to users: breaking update at version 0.5 (june 2013)
-------------
Up to version `<=0.5`, `stringdist` returned `-1` to indicate that either

* `maxDist` is exceeded, or
* the distance is undefined between two input strings

*From version `>=0.6`, this will be replaced by `Inf`.* The main reason:

* This allows easier comparison using `<` and `<=`

The old convention stemmed from the fact that some distances (e.g. q-gram, hamming) are strictly
integers. I'm abandoning this expressing all distances as `double` (R `numeric`) which allows me
to use `Inf`.


Installation
------------
To install the latest release from CRAN, open an R terminal and type

`install.packages('stringdist')`

To obtain the package from source code open a `bash` terminal (or `git bash` if you work under Windows
with `msysgit`) and type

```
git clone https://github.com/markvanderloo/stringdist.git
cd stringdist
bash ./build.bash
R CMD INSTALL output/stringdist_*.tar.gz
```

Warning: the github version can change any time and may not even build properly. As most
of the code is written in `C`, the development version may crash your `R`-session.


