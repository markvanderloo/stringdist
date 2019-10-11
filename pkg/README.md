
## stringdist

* Approximate matching and string distance calculations for R. 
* All distance and matching operations are system- and encoding-independent.
* Built for speed, using [openMP](https://www.openmp.org/) for parallel computing.

The package offers the following main functions:

* `stringdist`  computes pairwise distances between two input character vectors (shorter one is recycled)
* `stringdistmatrix` computes the distance matrix for one or two vectors
* `stringsim` computes a string similarity between 0 and 1, based on `stringdist`
* `amatch` is a fuzzy matching equivalent of R's native `match` function
* `ain` is a fuzzy matching equivalent of R's native `%in%` operator
* `seq_dist`, `seq_distmatrix`, `seq_amatch` and `seq_ain` for distances between, and matching of integer sequences. (see also the [hashr](https://github.com/markvanderloo/hashr) package).

These functions are built upon `C`-code that re-implements some common (weighted) string
distance functions. Distance functions include:

* Hamming distance; 
* Levenshtein distance (weighted)
* Restricted Damerau-Levenshtein distance (weighted, a.k.a. Optimal String Alignment)
* Full Damerau-Levenshtein distance
* Longest Common Substring distance
* Q-gram distance
* cosine distance for q-gram count vectors (= 1-cosine similarity)
* Jaccard distance for q-gram count vectors (= 1-Jaccard similarity)
* Jaro, and Jaro-Winkler distance
* Soundex-based string distance

Also, there are some utility functions:

* `qgrams()` tabulates the qgrams in one or more `character` vectors.
* `seq_qrams()` tabulates the qgrams (somtimes called ngrams) in one or more `integer` vectors.
* `phonetic()` computes phonetic codes of strings (currently only soundex)
* `printable_ascii()` is a utility function that detects non-printable ascii or non-ascii characters.

#### C API

As of version `0.9.5.0`  you can call a number of `stringdist` functions directly
from the `C` code of your R package. The description of the API can be found 

- By typing `?stringdist_api` in the R console
- Or open the vignette directly:

```
vignette("stringdist_C-Cpp_api", package="stringdist")
```

Examples of packages that link to `stringdist` can be found
[here](https://github.com/markvanderloo/linkstringdist) and
[here](https://github.com/ChrisMuir/refinr).


#### Resources

* A [paper](http://journal.r-project.org/archive/2014-1/loo.pdf) on stringdist has been published in the R-journal
* [Slides](http://www.slideshare.net/MarkVanDerLoo/stringdist-use-r2014) of te _useR!2014_ conference.

