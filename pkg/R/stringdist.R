#' A package for string distance calculation
#' 
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
#' @import parallel
{}


  
#' Compute distance metrics between strings
#'
#' @section Details:
#' \code{stringdist} computes pairwise string distances between elements of \code{character} vectors 
#' \code{a} and \code{b}, where the vector with less elements is recycled. 
#' 
#' \code{stringdistmatrix} computes the string distance matrix with rows according to
#' \code{a} and columns according to \code{b}.
#' 
#' Currently, the following distance metrics are supported:
#' \tabular{ll}{
#'    \code{osa} \tab Optimal string aligment, (restricted Damerau-Levenshtein distance).\cr
#'    \code{lv} \tab Levenshtein distance (as in R's native \code{\link[utils]{adist}}).\cr
#'    \code{dl} \tab Full Damerau-Levenshtein distance.\cr
#'    \code{hamming}  \tab Hamming distance (\code{a} and \code{b} must have same nr of characters).\cr
#'    \code{lcs} \tab Longest common substring distance.\cr
#'    \code{qgram} \tab \eqn{q}-gram distance. \cr
#'    \code{cosine} \tab cosine distance between \eqn{q}-gram profiles \cr
#'    \code{jaccard} \tab Jaccard distance between \eqn{q}-gram profiles \cr
#'    \code{jw} \tab Jaro, or Jaro-Winker distance.\cr
#'    \code{soundex} \tab Distance based on soundex encoding (see below)
#' }
#' 
#' Precise descriptions of the algorithms are given in the R-journal paper (see Citation section).
#' Below are some concise descriptions.
#' 
#' 
#' The \bold{Hamming distance} (\code{hamming}) counts the number of 
#' character substitutions that turns \code{b} into \code{a}. If \code{a} 
#' and \code{b} have different number of characters or if \code{maxDist} is 
#' exceeded, \code{Inf} is returned.
#'
#' The \bold{Levenshtein distance} (\code{lv}) counts the number of 
#' deletions, insertions and substitutions necessary to turn \code{b} into 
#' \code{a}. This method is equivalent to \code{R}'s native \code{\link[utils]{adist}} 
#' function. If \code{maxDist} is exceeded \code{Inf}  is returned.
#'
#' The \bold{Optimal String Alignment distance} (\code{osa}) is like the Levenshtein 
#' distance but also allows transposition of adjacent characters. Here, each 
#' substring  may be edited only once. (For example, a character cannot be transposed twice
#' to move it forward in the string). If \code{maxDist} is exceeded \code{Inf}  is returned.
#'
#' The \bold{full Damerau-Levensthein distance} (\code{dl}) allows for multiple 
#' edits on substrings. If \code{maxDist} is exceeded \code{Inf}  is returned.
#'
#' The \bold{longest common substring} is defined as the longest string that can be 
#' obtained by pairing characters from \code{a} and \code{b} while keeping the order 
#' of characters intact. The \bold{lcs-distance} is defined as the number of unpaired characters. 
#' The distance is equivalent to the edit distance allowing only deletions and insertions, 
#' each with weight one. If \code{maxDist} is exceeded \code{Inf}  is returned.
#'
#' A \bold{\eqn{q}-gram} is a subsequence of \eqn{q} \emph{consecutive} 
#' characters of a string. If \eqn{x} (\eqn{y}) is the vector of counts
#' of \eqn{q}-gram occurrences in \code{a} (\code{b}), the \bold{\eqn{q}-gram distance} 
#' is given by the sum over the absolute differences \eqn{|x_i-y_i|}.
#' The computation is aborted when \code{q} is is larger than the length of 
#' any of the strings. In that case \code{Inf}  is returned.
#'
#' The \bold{cosine distance} is computed as \eqn{1-x\cdot y/(\|x\|\|y\|)}, where \eqn{x} and 
#' \eqn{y} were defined above.
#' 
#' Let \eqn{X} be the set of unique \eqn{q}-grams in \code{a} and \eqn{Y} the set of unique 
#' \eqn{q}-grams in \code{b}. The \bold{Jaccard distance} is given by \eqn{1-|X\cap Y|/|X\cup Y|}.
#'
#' The \bold{Jaro distance} (\code{method='jw'}, \code{p=0}), is a number
#' between 0 (exact match) and 1 (completely dissimilar) measuring 
#' dissimilarity between strings.  It is defined to be 0 when both strings have
#' length 0, and 1 when  there are no character matches between \code{a} and
#' \code{b}.  Otherwise, the Jaro distance is defined as 
#' \eqn{1-(1/3)(w_1m/|a| + w_2m/|b| + w_3(m-t)/m)}. 
#' Here,\eqn{|a|} indicates the number of characters in \code{a}, \eqn{m} is
#' the number of character matches and \eqn{t} the number of transpositions of
#' matching characters. The \eqn{w_i} are weights associated with the characters
#' in \code{a}, characters in \code{b} and with transpositions.  A character
#' \eqn{c} of \code{a} \emph{matches} a character from \code{b} when \eqn{c}
#' occurs in \code{b}, and the index of \eqn{c} in \code{a} differs less than
#' \eqn{\max(|a|,|b|)/2 -1} (where we use integer division) from the index of
#' \eqn{c} in \code{b}. Two matching characters are transposed when they are
#' matched but they occur in different order in string \code{a} and \code{b}.
#'  
#' The \bold{Jaro-Winkler distance} (\code{method=jw}, \code{0<p<=0.25}) adds a
#' correction term to the Jaro-distance. It is defined as \eqn{d - l*p*d}, where
#' \eqn{d} is the Jaro-distance. Here,  \eqn{l} is obtained by counting, from
#' the start of the input strings, after how many characters the first
#' character mismatch between the two strings occurs, with a maximum of four. The
#' factor \eqn{p} is a penalty factor, which in the work of Winkler is often
#' chosen \eqn{0.1}.
#'
#' For the \bold{soundex} method, strings are translated to a soundex code (see \code{\link{phonetic}} for a specification). The
#' distance between strings is 0 when they have the same soundex code,
#' otherwise 1. Note that soundex recoding is only meaningful for characters
#' in the ranges a-z and A-Z. A warning is emitted when non-printable or non-ascii
#' characters are encountered. Also see \code{\link{printable_ascii}}.
#'
#' @section Encoding issues:
#' If \code{bytes=FALSE}, input strings are re-encoded to \code{utf8} an then to \code{integer}
#' vectors prior to the distance calculation (since the underlying \code{C}-code expects \code{unsigned int}s). 
#' This double conversion is necessary as it seems the only way to
#' reliably convert (possibly multibyte) characters to integers on all systems
#' supported by \code{R}. \code{R}'s native \code{\link[utils]{adist}} function does this as well. 
#'
#' If \code{bytes=TRUE}, the input strings are treated as if each byte was a
#' single character. This may be significantly faster since it avoids
#' conversion of \code{utf8} to integer with \code{\link[base]{utf8ToInt}} (up to a factor of 3, for strings of 5-25 characters).
#' However, results may depend on the (possibly multibyte)  character encoding scheme
#' and note that \code{R}'s internal encoding scheme is OS-dependent. 
#' If you're sure that all your input is \code{ASCII},  you can safely set 
#' \code{useBytes=TRUE} to profit from the speed gain on any platform. 
#'
#' See base \code{R}'s \code{\link[base]{Encoding}} and
#' \code{\link[base]{iconv}} documentation for details on how \code{R} handles
#' character encoding. 
#'
#' @section Unicode normalisation:
#' In \code{utf-8}, the same (accented) character may be represented as several byte sequences. For example, an u-umlaut
#' can be represented with a single byte code or as a byte code representing \code{'u'} followed by a modifier byte code
#' that adds the umlaut. The \href{http://cran.r-project.org/web/packages/stringi/}{stringi} package 
#' of Gagolevski and Tartanus offers unicode normalisation tools. 
#'
#'
#' @section Paralellization:
#' The \code{stringdistmatrix} function uses \code{\link[parallel]{makeCluster}} to create a local cluster and compute the
#' distance matrix in parallel when \code{ncores>1}. The cluster is terminated after the matrix has been computed. 
#' As the cluster is local, the \code{ncores} parameter should not be larger than the number
#' of cores on your machine. Use \code{\link[parallel]{detectCores}} to check the number of cores available. Alternatively,
#' you can create a cluster using \code{\link[parallel]{makeCluster}} and pass that to \code{stringdistmatrix} (through the \code{cluster} argument. 
#' This allows you to reuse the cluster setup for other calculations.
#' There is overhead in creating clusters, so creating the cluster yourself is a good choice if you want to call \code{stringdistmatrix} 
#' multiple times, for example in a loop.
#'
#' @section Citation:
#' If you would like to cite this package, please cite the R-journal paper: 
#' \itemize{
#' \item{M.P.J. van der Loo (2014). The \code{stringdist} package for approximate string matching. 
#'  R Journal 6(1) pp 111-122}
#' }
#' Or use \code{citation('stringdist')} to get a bibtex item.
#'
#' @references
#'
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
#'  An extensive overview of offline string matching algorithms is given by L. Boytsov (2011). Indexing
#'  methods for approximate dictionary searching: comparative analyses. ACM Journal of experimental
#'  algorithmics 16 1-88.
#' }
#' \item{
#'  An extensive overview of (online) string matching algorithms is given by G. Navarro (2001). 
#'  A guided tour to approximate string matching, ACM Computing Surveys 33 31-88.
#' }
#' \item{
#' Many algorithms are available in pseudocode from wikipedia: \url{http://en.wikipedia.org/wiki/Damerau-Levenshtein_distance}.
#' }
#' \item{The code for the full Damerau-Levenshtein distance was adapted from Nick Logan's
#'  \href{https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c}{public github repository}.
#' }
#'
#' \item{
#' A good reference for qgram distances is E. Ukkonen (1992), Approximate string matching with q-grams and maximal matches. 
#' Theoretical Computer Science, 92, 191-211.
#' }
#'
#' \item{\href{http://en.wikipedia.org/wiki/Jaro\%E2\%80\%93Winkler_distance}{Wikipedia} describes the Jaro-Winker
#' distance used in this package. Unfortunately, there seems to be no single
#'  definition for the Jaro distance in literature. For example Cohen, Ravikumar and Fienberg (Proceeedings of IIWEB03, Vol 47, 2003)
#'  report a different matching window for characters in strings \code{a} and \code{b}. 
#' }
#'
#' \item{Raffael Vogler wrote a nice 
#' \href{http://www.joyofdata.de/blog/comparison-of-string-distance-algorithms/}{blog}
#' comparing different string distances in this package.
#'
#'
#'}
#'
#'}
#'
#'
#'
#' @param a R object (target); will be converted by \code{as.character}.
#' @param b R object (source); will be converted by \code{as.character}.
#' @param method Method for distance calculation. The default is \code{"osa"} (see details).
#' @param useBytes Perform byte-wise comparison. \code{useBytes=TRUE} is faster but may yield different
#' 	results depending on character encoding. See also below, under ``encoding issues''.
#' @param weight For \code{method='osa'} or \code{'dl'}, the penalty for deletion, insertion, substitution and transposition, in that order.
#'   When \code{method='lv'}, the penalty for transposition is ignored. When \code{method='jw'}, the weights associated with characters
#'   of \code{a}, characters from \code{b} and the transposition weight, in that order.
#'   Weights must be positive and not exceed 1. \code{weight} is
#'   ignored completely when \code{method='hamming'}, \code{'qgram'}, \code{'cosine'}, \code{'Jaccard'}, or \code{'lcs'}. 
#' @param maxDist  [DEPRECATED AND WILL BE REMOVED] Currently kept for backward compatibility. It does not
#' offer any speed gain. (In fact, it currently slows things down when set to anything different from \code{Inf}).
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param p Penalty factor for Jaro-Winkler distance. The valid range for \code{p} is \code{0 <= p <= 0.25}. 
#'  If \code{p=0} (default), the Jaro-distance is returned. Applies only to \code{method='jw'}.
#' @param nthread (positive integer) Number of threads to use.
#'
#'
#' @return For \code{stringdist},  a vector with string distances of size \code{max(length(a),length(b))}.
#'  For \code{stringdistmatrix}, a \code{length(a)xlength(b)} \code{matrix}. The returned distance is
#'  nonnegative if it can be computed, \code{NA} if any of the two argument strings is \code{NA} and \code{Inf}
#'  when it cannot be computed or \code{maxDist} is exceeded. See details for the meaning of \code{Inf} for the various algorithms.
#'  
#'  
#' @example ../examples/stringdist.R
#' @export
stringdist <- function(a, b
  , method=c("osa","lv","dl","hamming","lcs", "qgram","cosine","jaccard","jw","soundex")
  , useBytes = FALSE
  , weight=c(d=1,i=1,s=1,t=1) 
  , maxDist=Inf, q=1, p=0
  , nthread = 1L
){

  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 || length(b) == 0){ 
    return(numeric(0))
  }
  if ( max(length(a),length(b)) %% min(length(a),length(b)) != 0 ){
    warning(RECYCLEWARNING)
  }
  method <- match.arg(method)
  nthread <- as.integer(nthread)
  stopifnot(
      all(is.finite(weight))
      , all(weight > 0)
      , all(weight <=1)
      , q >= 0
      , p <= 0.25
      , p >= 0
      , is.logical(useBytes)
      , ifelse(method %in% c('osa','dl'), length(weight) >= 4, TRUE)
      , ifelse(method %in% c('lv','jw') , length(weight) >= 3, TRUE)
      , nthread > 0
  )
  if (!useBytes){
    a <- char2int(a)
    b <- char2int(b)
  }
  if (method == 'jw') weight <- weight[c(2,1,3)]
  do_dist(b, a, method, weight, maxDist, q, p, nthread)
}


#' @param useNames Use input vectors as row and column names?
#' @param ncores Number of cores to use. If \code{ncores>1}, a local cluster is
#' created using \code{\link[parallel]{makeCluster}}. Parallelisation is over \code{b}, so 
#' the speed gain by parallelisation is highest when \code{b} has less elements than \code{a}.
#' @param cluster (Optional) a custom cluster, created with
#' \code{\link[parallel]{makeCluster}}. If \code{cluster} is not \code{NULL},
#' \code{ncores} is ignored.
#'
#'
#'
#' @rdname stringdist
#' @export
stringdistmatrix <- function(a, b, 
  method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard","jw","soundex"), 
  useBytes = FALSE,
  weight=c(d=1,i=1,s=1,t=1),  maxDist=Inf, q=1, p=0,
  useNames=FALSE, ncores=1, cluster=NULL
){
  
  a <- as.character(a)
  b <- as.character(b)
 
  if (length(a) == 0 || length(b) == 0){ 
   return(matrix(numeric(0)))
  }

  if (useNames){
   rowns <- a
   colns <- b
  }

  method <- match.arg(method)
  stopifnot(
      all(is.finite(weight))
      , all(weight > 0)
      , all(weight <=1)
      , q >= 0
      , p <= 0.25
      , p >= 0
      , is.logical(useBytes)
      , ifelse(method %in% c('osa','dl'), length(weight) >= 4, TRUE)
      , ifelse(method %in% c('lv','jw') , length(weight) >= 3, TRUE)
      , ncores > 0
  )
  if (!useBytes){
    a <- char2int(a)
    b <- lapply(char2int(b),list)
  }
  if (method == 'jw') weight <- weight[c(2,1,3)]
  if (ncores==1){
    x <- sapply(b,do_dist, USE.NAMES=FALSE, a,method,weight,maxDist, q, p)
  } else {
    if ( is.null(cluster) ){
      cluster <- makeCluster(ncores)
      turn_cluster_off <- TRUE
    } else {
      stopifnot(inherits(cluster, 'cluster'))
      turn_cluster_off <- FALSE
    }
    x <- parSapply(cluster, b,do_dist,a,method,weight,maxDist, q, p)
    if (turn_cluster_off) stopCluster(cluster)
  }
  if (!useNames)  as.matrix(x) else structure(as.matrix(x), dimnames=list(rowns,colns))
}


char2int <- function(x){
  # For some OS's enc2utf8 had unexpected behavior for NA's,
  # see https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15201.
  # This is fixed for R >= 2.15.3.
  # i <- !is.na(x)
  # x[i] <- enc2utf8(x[i])
  lapply(enc2utf8(x),utf8ToInt)
}



do_dist <- function(a, b, method, weight, maxDist, q, p, nthread=1L){
  d <- switch(method,
    osa     = .Call('R_osa'   , a, b, as.double(weight)),
    lv      = .Call('R_lv'    , a, b, as.double(weight)),
    dl      = .Call('R_dl'    , a, b, as.double(weight)),
    hamming = .Call('R_hm'    , a, b, nthread),
    lcs     = .Call('R_lcs'   , a, b),
    qgram   = .Call('R_qgram_tree' , a, b, as.integer(q), 0L),
    cosine  = .Call('R_qgram_tree' , a, b, as.integer(q), 1L),
    jaccard = .Call('R_qgram_tree' , a, b, as.integer(q), 2L),
    jw      = .Call('R_jw'    , a, b, as.double(p), as.double(weight)),
    soundex = .Call('R_soundex_dist', a, b)
  )
  if (maxDist < Inf ){
    d[!is.na(d) & d > maxDist] <- Inf
  }
  d
}



