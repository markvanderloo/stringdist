#' A package for string distance calculation
#' 
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
{}

#' Compute distance between strings
#'
#' @section Details:
#' \code{stringdist} computes pairwise string distances between elements of \code{character} vectors \code{a} and \code{b},
#' where the vector with less elements is recycled. \code{stringdistmatrix} computes the string distance matrix with rows according to
#' \code{a} and columns according to \code{b}.
#'
#' 
#' Currently, the following distance metrics are supported:
#' \tabular{ll}{
#'    \code{osa} \tab Optimal string aligment, (restricted Damerau-Levenshtein distance).\cr
#'    \code{lv} \tab Levenshtein distance.\cr
#'    \code{dl} \tab Full Damerau-Levenshtein distance.\cr
#'    \code{h}  \tab Hamming distance (\code{a} and \code{b} must have same nr of characters).\cr
#'    \code{lcs} \tab Longest common substring.\cr
#'    \code{qgram} \tab \eqn{q}-gram distance. \cr
#'    \code{jw} \tab Jaro-Winker distance.
#' }
#' The \bold{Hamming distance} counts the number of character substitutions that turns 
#' \code{b} into \code{a}. If \code{a} and \code{b} have different number of characters \code{-1} is
#' returned.
#'
#' The \bold{Levenshtein distance} (\code{lv}) counts the number of deletions, insertions and substitutions necessary
#' to turn \code{b} into \code{a}. This method is equivalent to \code{R}'s native \code{adist} function.
#' The computation is aborted when \code{maxDist} is exceeded, in which case \code{-1}  is returned.
#'
#' The \bold{Optimal String Alignment distance} (\code{osa}) is like the Levenshtein distance but also 
#' allows transposition of adjacent characters. Here, each substring  may be edited only once so a 
#' character cannot be transposed twice. 
#' The computation is aborted when \code{maxDist} is exceeded, in which case \code{-1}  is returned.
#'
#' The \bold{full Damerau-Levensthein distance} (\code{dl}) allows for multiple transpositions.
#' The computation is aborted when \code{maxDist} is exceeded, in which case \code{-1}  is returned.
#'
#' The \bold{longest common substring} is defined as the longest string that can be obtained by pairing characters
#' from \code{a} and \code{b} while keeping the order of characters intact. The lcs-distance is defined as the
#' number of unpaired characters. The distance is equivalent to the edit distance allowing only deletions and
#' insertions, each with weight one.
#' The computation is aborted when \code{maxDist} is exceeded, in which case \code{-1}  is returned.
#'
#' A \bold{\eqn{q}-gram} is a subsequence of \eqn{q} \emph{consecutive} characters of a string. If \eqn{x} (\eqn{y}) is the vector of counts
#' of \eqn{q}-gram occurrences in \code{a} (\code{b}), the \eqn{q}-gram distance is given by the sum over
#' the absolute differences \eqn{|x_i-y_i|}.
#' The computation is aborted when \code{q} is is larger than the length of any of the strings. In that case \code{-1}  is returned.
#'
#' The \bold{Jaro distance} (\code{method=jw}, \code{p=0}), is a number between 0 (exact match) and 1 (completely dissimilar) measuring 
#' dissimilarity between strings.
#' It is defined to be 0 when both strings have length 0, and 1 when  there are no character matches between \code{a} and \code{b}. 
#' Otherwise, the Jaro distance is defined as \eqn{1-(1/3)(m/|a| + m/|b| + (m-t)/m)}. Here,\eqn{|a|} indicates the number of
#' characters in \code{a} (after conversion to integers), \eqn{m} is the number of 
#' character matches and \eqn{t} the number of transpositions of matching characters.
#' A character \eqn{c} of \code{a} \emph{matches} a character from \code{b} when
#' \eqn{c} occurs in \code{b}, and the index of \eqn{c} in \code{a} differs less than \eqn{max(|a|,|b|)/2 -1} (where we use integer division).
#' Two matching characters are transposed when they are matched but they occur in different order in string \code{a} and \code{b}.
#'  
#' The \bold{Jaro-Winkler distance} (\code{method=jw}, \code{0<p<=0.25}) adds a correction term to the Jaro-distance. It is defined as \eqn{d - l*p*d}, where
#' \eqn{d} is the Jaro-distance. Here,  \eqn{l} is obtained by counting, from the start of the input strings, after how many
#' characters the first character mismatch between the two strings occurs, with a maximum of four. The factor \eqn{p}
#' is a penalty factor, which in the work of Winkler is often chosen \eqn{0.1}.
#'
#' @section Encoding issues:
#' Input strings are re-encoded to \code{utf8} an then to \code{integer}
#' vectors prior to the distance calculation (since the underlying \code{C}-code expects unsigned ints). 
#' This double conversion is necessary as it seems the only way to
#' reliably convert (possibly multibyte) characters to integers on all systems
#' supported by \code{R}.
#' (\code{R}'s native \code{\link[utils]{adist}} function does this as well). 
#' See \code{\link[base]{Encoding}} for further details.
#'
#' @section Paralellization:
#' The \code{stringdistmatrix} function uses \code{\link[parallel]{makeCluster}} to generate a cluster and compute the
#' distance matrix in parallel.  As the cluster is local, the \code{ncores} parameter should not be larger than the number
#' of cores on your machine. Use \code{\link[parallel]{detectCores}} to check the number of cores available. Alternatively,
#' you can create a cluster by yourself, using \code{\link[parallel]{makeCluster}} and pass that to \code{stringdistmatrix}.
#'
#' @references
#' \itemize{
#' \item{
#'    R.W. Hamming (1950). Error detecting and Error Correcting codes, The Bell System Technical Journal 29, 147-160
#'  }
#' \item{
#'  V.I. Levenshtein. (1960). Binary codes capable of correcting deletions, insertions, and reversals. Soviet Physics Doklady 10 707-711.
#' }
#' \item{
#' F.J. Damerau (1964) A technique for computer detection and correction of spelling errors. Communications of the ACM 7 171-176.
#' }
#' \item{
#' Many algorithms are available in pseudocode from wikipedia: http://en.wikipedia.org/wiki/Damerau-Levenshtein_distance.
#' }
#' \item{The code for the full Damerau-Levenshtein distance was adapted from Nick Logan's public github repository:
#'  \url{https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c}.
#' }
#'
#' \item{
#' A good reference for qgram distances is E. Ukkonen (1992), Approximate string matching with q-grams and maximal matches. 
#' Theoretical Computer Science, 92, 191-211.
#' }
#'
#' \item{Wikipedia \code{http://en.wikipedia.org/wiki/Jaro\%E2\%80\%93Winkler_distance} describes the Jaro-Winker
#' distance used in this package. Unfortunately, there seems to be no single
#'  definition for the Jaro distance in literature. For example Cohen, Ravikumar and Fienberg (Proceeedings of IIWEB03, Vol 47, 2003)
#'  report a different matching window for characters in strings \code{a} and \code{b}. 
#' }
#'
#'}
#'
#'
#'
#' @param a R object (target); will be converted by \code{as.character}.
#' @param b R object (source); will be converted by \code{as.character}.
#' @param method Method for distance calculation. The default is \code{"osa"} (see details).
#' @param weight The penalty for deletion, insertion, substitution and transposition, in that order.  
#'   Weights must be positive and not exceed 1. \code{weight[4]} is ignored when \code{method='lv'} and \code{weight} is
#'   ignored completely when \code{method='h'}, \code{method='qgram'} or \code{method='lcs'}.
#' @param maxDist  Maximum string distance before calculation is stopped, \code{maxDist=0} 
#'    means calculation goes on untill the distance is computed. Ignored for \code{method='qgram'} and
#'    \code{method='jw'}.
#' @param q  size of the \eqn{q}-gram, must be nonnegative. Ignored for all but \code{method='qgram'}.
#' @param p penalty factor for Jaro-Winkler distance. The valid range for \code{p} is \code{0<= p <= 0.25}. 
#'  If \code{p=0} (default), the Jaro-distance is returned. Ignored for all methods except \code{'jw'}.
#'
#'
#'
#'
#'
#'
#' @return For \code{stringdist},  a vector with string distances of size \code{max(length(a),length(b))}.
#'  For \code{stringdistmatrix}, a \code{length(a)xlength(b)} \code{matrix}. The returned distance is
#'  nonnegative if it can be computed, \code{NA} if any of the two argument strings is \code{NA} and \code{-1}
#'  when it cannot be computed. See details for the meaning of \code{-1} for the various algorithms.
#'  
#'  
#' @example ../examples/stringdist.R
#' @export
stringdist <- function(a, b, 
  method=c("osa","lv","dl","h","lcs", "qgram", "jw"), 
  weight=c(d=1,i=1,s=1,t=1), 
  maxDist=0, q=1, p=0
){
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 || length(b) == 0){ 
    return(numeric(0))
  }
  method <- match.arg(method)
  a <- char2int(a)
  b <- char2int(b)
  stopifnot(
      all(is.finite(weight)),
      all(weight > 0),
      all(weight <=1),
      p <= 0.25,
      p >= 0
  )
  do_dist(b,a,method,weight,maxDist,q,p)
}


#' @param ncores number of cores to use. If \code{ncores>1}, a local cluster is
#' created using \code{\link[parallel]{makeCluster}}. Parallelisation is over \code{b}, so 
#' the speed gain by parallelisation is highest when \code{b} has less elements than \code{a}.
#' @param cluster (optional) a custom cluster, created with
#' \code{\link[parallel]{makeCluster}}. If \code{cluster} is not \code{NULL},
#' \code{ncores} is ignored.
#'
#'
#'
#' @rdname stringdist
#' @export
stringdistmatrix <- function(a, b, 
  method=c("osa","lv","dl","h","lcs","qgram", "jw"), 
  weight=c(d=1,i=1,s=1,t=1), 
  maxDist=0, q=1, p=0,
  ncores=1, cluster=NULL
){
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 || length(b) == 0){ 
   return(numeric(0))
  }
  method <- match.arg(method)
  stopifnot(
      all(is.finite(weight)),
      all(weight > 0),
      all(weight <=1),
      p <= 0.25,
      p >= 0
  )
  a <- char2int(a)
  b <- lapply(char2int(b),list)
  if (ncores==1){
    x <- sapply(b,do_dist,a,method,weight,maxDist, q)
  } else {
    if ( is.null(cluster) ){
      cl <- makeCluster(ncores)
    } else {
      stopifnot(inherits(cluster, 'cluster'))
      cl <- cluster
    }
    x <- parSapply(cluster, b,do_dist,a,method,weight,maxDist, q, p)
    if (is.null(cluster)) stopCluster(cl)
  }
  x
}


char2int <- function(x){
  # For some OS's enc2utf8 has unexpected behavior for NA's,
  # see https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15201.
  # This is fixed for R >= 2.15.3.
  # i <- !is.na(x)
  # x[i] <- enc2utf8(x[i])
  lapply(enc2utf8(x),utf8ToInt)
}



do_dist <- function(a, b, method, weight, maxDist, q, p){
  switch(method,
    osa     = .Call('R_osa'   , a, b, as.double(weight), as.double(maxDist)),
    lv      = .Call('R_lv'    , a, b, as.double(weight), as.double(maxDist)),
    dl      = .Call('R_dl'    , a, b, as.double(weight), as.double(maxDist)),
    h       = .Call('R_hm'    , a, b, as.integer(maxDist)),
    lcs     = .Call('R_lcs'   , a, b, as.integer(maxDist)),
    qgram   = .Call('R_qgram_tree' , a, b, as.integer(q)),
    jw      = .Call('R_jaro_winkler', a, b, as.double(p))
  )
}



