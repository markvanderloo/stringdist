#' A package for string distance calculation and approximate string matching.
#'
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
#' @import parallel
#'
#' @section Supported distances:
#'
#'
#' @section Citation:
#' If you would like to cite this package, please cite the \href{http://journal.r-project.org/archive/2014-1/loo.pdf}{R Journal Paper}: 
#' \itemize{
#' \item{M.P.J. van der Loo (2014). The \code{stringdist} package for approximate string matching. 
#'  R Journal 6(1) pp 111-122}
#' }
#' Or use \code{citation('stringdist')} to get a bibtex item.
#'
#'
#' 
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
#' 
#' A short description of these algorithms is provided \link[=stringdist-package]{here}, or
#' see the \href{http://journal.r-project.org/archive/2014-1/loo.pdf}{R Journal Paper} (external link) for
#' more formal descriptions.
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
#'
#' By default \code{stringdist} and \code{stringdistmatrix} will use \code{getOption("sd_num_thread")} threads.
#' When the package is loaded, this option is set to the smaller of the number of available cores or the
#' environment variable \code{OMP_THREAD_LIMIT}, if available. The number of cores is detected with
#' \code{parallel::detectCores}. Using the maximum number of threads is not allways the fastest option.
#' At least one core will also be occupied with for example OS services, so it may be faster to use one core less
#' than the maximum number of cores.
#' 
#' In older versions (<0.9) of \code{stringdist}, the \code{cluster} and \code{ncores} argument were the only 
#' paralellization options, and only for \code{stringdistmatrix}. These options are based on the parallel package 
#' which starts multiple R-sessions to run R code in parallel. If you're running R on a single machine it is both 
#' faster and easier to use the default multithreading (so do not specify \code{ncores} or \code{cluster}).
#'
#' As of the introduction of the \code{nthreads} argument, the \code{ncores} is mostly useless, although it still works.
#' If \code{ncores>0}, a local cluster of running R-sessions is set up automatically. Each R-session will use \code{nthread}
#' threads.
#'
#' The \code{cluster} argument is only interesting when the cluster is set up over different physical nodes. For example when
#' setting up a network of nodes accross physically different machines. In each node, \code{nthread} threads will be used.
#'
#'
#' @section Acknowledgements:
#' The code for the full Damerau-Levenshtein distance was adapted from Nick Logan's
#' \href{https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c}{public github repository}.
#' 
#' C code for converting UTF-8 to integer was copied from the R core.
#'
#' @section Other:
#'
#' \itemize{
#'
#' \item{Raffael Vogler wrote a nice 
#' \href{http://www.joyofdata.de/blog/comparison-of-string-distance-algorithms/}{blog}
#' comparing different string distances in this package.
#'}
#'
#'}
#'
#'
#'
#' @param a R object (target); will be converted by \code{as.character}.
#' @param b R object (source); will be converted by \code{as.character}.
#' @param method Method for distance calculation. The default is \code{"osa"}, see \code{\link{stringdist-metrics}}.
#' @param useBytes Perform byte-wise comparison. \code{useBytes=TRUE} is faster but may yield different
#' 	results depending on character encoding. See also below, under ``encoding issues''.
#' @param weight For \code{method='osa'} or \code{'dl'}, the penalty for deletion, insertion, substitution and transposition, in that order.
#'   When \code{method='lv'}, the penalty for transposition is ignored. When \code{method='jw'}, the weights associated with characters
#'   of \code{a}, characters from \code{b} and the transposition weight, in that order.
#'   Weights must be positive and not exceed 1. \code{weight} is
#'   ignored completely when \code{method='hamming'}, \code{'qgram'}, \code{'cosine'}, \code{'Jaccard'}, \code{'lcs'}, or \code{soundex}. 
#' @param maxDist  [DEPRECATED AND WILL BE REMOVED] Currently kept for backward compatibility. It does not
#' offer any speed gain. (In fact, it currently slows things down when set to anything different from \code{Inf}).
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param p Penalty factor for Jaro-Winkler distance. The valid range for \code{p} is \code{0 <= p <= 0.25}. 
#'  If \code{p=0} (default), the Jaro-distance is returned. Applies only to \code{method='jw'}.
#' @param nthread (positive integer) Number of threads used by the underlying C-code. The default is the number of cores
#'  detected by \code{\link[parallel]{detectCores}}.
#'
#'
#' @return For \code{stringdist},  a vector with string distances of size \code{max(length(a),length(b))}.
#'  For \code{stringdistmatrix}, a \code{length(a)xlength(b)} \code{matrix}. The returned distance is
#'  nonnegative if it can be computed, \code{NA} if any of the two argument strings is \code{NA} and \code{Inf}
#'  when \code{maxDist} is exceeded or, in case of the hamming distance, when the two compared strings have different length.
#'  
#'  
#' @example ../examples/stringdist.R
#' @export
stringdist <- function(a, b
  , method=c("osa","lv","dl","hamming","lcs", "qgram","cosine","jaccard","jw","soundex")
  , useBytes = FALSE
  , weight=c(d=1,i=1,s=1,t=1) 
  , maxDist=Inf, q=1, p=0
  , nthread = getOption("sd_num_thread")
){
  if (maxDist < Inf)
    message("Argument 'maxDist' is deprecated for function 'stringdist'")
  # note: enc2utf8 is very efficient when the native encoding is already UTF-8.
  a <- as.character(a)
  b <- as.character(b)
  if ( !useBytes ){
    a <- enc2utf8(a)
    b <- enc2utf8(b)
  }

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
  if (method == 'jw') weight <- weight[c(2,1,3)]
  do_dist(b, a, method, weight, maxDist, q, p, useBytes, nthread)
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
stringdistmatrix <- function(a, b
  , method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard","jw","soundex")
  , useBytes = FALSE
  , weight=c(d=1,i=1,s=1,t=1),  maxDist=Inf, q=1, p=0
  , useNames=FALSE, ncores=1, cluster=NULL
  , nthread = getOption("sd_num_thread")
){
  if (maxDist < Inf)
    message("Argument 'maxDist' is deprecated for function 'stringdistmatrix'") 
 
  a <- as.character(a)
  b <- as.character(b)

  if (!useBytes){
    a <- enc2utf8(a)
    b <- enc2utf8(b)
  }
 
  if (length(a) == 0 || length(b) == 0){ 
   return(matrix(numeric(0)))
  }

  if (useNames){
   rowns <- a
   colns <- b
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
      , ncores > 0
      , nthread > 0
  )
  if (method == 'jw') weight <- weight[c(2,1,3)]
  if (ncores==1){
    x <- sapply(b,do_dist, USE.NAMES=FALSE, a,method,weight,maxDist, q, p,useBytes, nthread)
  } else {
    if ( is.null(cluster) ){
      cluster <- makeCluster(ncores)
      turn_cluster_off <- TRUE
    } else {
      stopifnot(inherits(cluster, 'cluster'))
      turn_cluster_off <- FALSE
    }
    x <- parSapply(cluster, b,do_dist,a,method,weight,maxDist, q, p, useBytes, nthread)
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



do_dist <- function(a, b, method, weight, maxDist, q, p, useBytes=FALSE, nthread=1L){
  d <- switch(method,
    osa     = .Call('R_osa'   , a, b, as.double(weight), useBytes, nthread),
    lv      = .Call('R_lv'    , a, b, as.double(weight), useBytes, nthread),
    dl      = .Call('R_dl'    , a, b, as.double(weight), useBytes, nthread),
    hamming = .Call('R_hm'    , a, b, useBytes, nthread),
    lcs     = .Call('R_lcs'   , a, b, useBytes, nthread),
    qgram   = .Call('R_qgram_tree' , a, b, as.integer(q), 0L, useBytes, nthread),
    cosine  = .Call('R_qgram_tree' , a, b, as.integer(q), 1L, useBytes, nthread),
    jaccard = .Call('R_qgram_tree' , a, b, as.integer(q), 2L, useBytes, nthread),
    jw      = .Call('R_jw'    , a, b, as.double(p), as.double(weight), useBytes, nthread),
    soundex = .Call('R_soundex_dist', a, b, useBytes, nthread)
  )
  if (maxDist < Inf ){
    d[!is.na(d) & d > maxDist] <- Inf
  }
  d
}



