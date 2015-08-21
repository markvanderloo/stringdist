#' A package for string distance calculation and approximate string matching.
#'
#' @section Introduction:
#'
#' The \pkg{stringdist} package offers fast and platform-independent string
#' metrics. Its main purpose is to compute various string distances and to do 
#' approximate text matching between character vectors. As of version 0.9.3,
#' it is also possible to compute distances between sequences represented by
#' integer vectors.
#' 
#' 
#' A typical use is to match strings that are not precisely the same. For
#' example
#'
#' \code{  amatch(c("hello","g'day"),c("hi","hallo","ola"),maxDist=2)}
#'
#' returns \code{c(2,NA)} since \code{"hello"} matches closest with
#' \code{"hallo"}, and within the maximum (optimal string alignment) distance.
#' The second element, \code{"g'day"}, matches closest with \code{"ola"} but
#' since the distance equals 4, no match is reported.
#'
#' A second typical use is to compute string distances. For example 
#'
#' \code{  stringdist(c("g'day"),c("hi","hallo","ola"))}
#'
#' Returns \code{c(5,5,4)} since these are the distances between \code{"g'day"}
#' and respectively \code{"hi"}, \code{"hallo"}, and \code{"ola"}.
#'
#' A third typical use would be to compute a \code{dist} object. The command
#' 
#' \code{stringdistmatrix(c("foo","bar","boo","baz"))}
#'
#' returns an object of class \code{dist} that can be used by clustering
#' algorithms such as \code{stats::hclust}.
#' 
#' A fourth use is to compute string distances between general sequences,
#' represented as integer vectors (which must be stored in a \code{list}):
#'
#' \code{seq_dist( list(c(1L,1L,2L)), list(c(1L,2L,1L),c(2L,3L,1L,2L)) )}
#'
#' The above code yields the vector \code{c(1,2)} (the first shorter first
#' argument is recycled over the longer second argument)
#' 
#' Besides documentation for each function, the main topics documented are:
#' 
#' \itemize{
#' \item{\code{\link{stringdist-metrics}} -- string metrics supported by the package}
#' \item{\code{\link{stringdist-encoding}} -- how encoding is handled by the package}
#' \item{\code{\link{stringdist-parallelization}} -- on multithreading }
#' }
#'
#' @section Acknowledgements:
#' \itemize{
#'   \item{The code for the full Damerau-Levenshtein distance was adapted from Nick Logan's
#'   \href{https://github.com/ugexe/Text--Levenshtein--Damerau--XS/blob/master/damerau-int.c}{public github repository}.}
#'   \item{C code for converting UTF-8 to integer was copied from the R core for performance reasons.}
#'   \item{The code for soundex conversion and string similarity was kindly contributed by Jan van der Laan.}
#' }
#' @section Citation:
#' If you would like to cite this package, please cite the \href{http://journal.r-project.org/archive/2014-1/loo.pdf}{R Journal Paper}: 
#' \itemize{
#' \item{M.P.J. van der Loo (2014). The \code{stringdist} package for approximate string matching. 
#'  R Journal 6(1) pp 111-122}
#' }
#' Or use \code{citation('stringdist')} to get a bibtex item.
#'
#' @name stringdist-package
#' @docType package
#' @useDynLib stringdist
#' @importFrom parallel detectCores 
#'
#'
#' 
{}

listwarning <- function(x,y){
sprintf("
You are passing one or more arguments of type 'list' to
'%s'. These arguments will be converted with 'as.character'
which is likeley not to give what you want (did you mean to use '%s'?).
This warning can be avoided by explicitly converting the argument(s).
",x,y)
}

#' Compute distance metrics between strings
#'
#' 
#' \code{stringdist} computes pairwise string distances between elements of
#' \code{a} and \code{b}, where the argument with less elements is recycled.
#' \code{stringdistmatrix} computes the string distance matrix with rows
#' according to
#' \code{a} and columns according to \code{b}.
#' 
#'
#' @param a R object (target); will be converted by \code{as.character}
#' @param b R object (source); will be converted by \code{as.character}
#'   This argument is optional for \code{stringdistmatrix} (see section \code{Value}).
#' @param method Method for distance calculation. The default is \code{"osa"},
#'   see \code{\link{stringdist-metrics}}.
#' @param useBytes Perform byte-wise comparison, see
#'   \code{\link{stringdist-encoding}}.
#' @param weight For \code{method='osa'} or \code{'dl'}, the penalty for
#'   deletion, insertion, substitution and transposition, in that order. When
#'   \code{method='lv'}, the penalty for transposition is ignored. When
#'   \code{method='jw'}, the weights associated with characters of \code{a},
#'   characters from \code{b} and the transposition weight, in that order. 
#'   Weights must be positive and not exceed 1. \code{weight} is ignored
#'   completely when \code{method='hamming'}, \code{'qgram'}, \code{'cosine'},
#'   \code{'Jaccard'}, \code{'lcs'}, or \code{soundex}.
#' @param maxDist  [DEPRECATED AND WILL BE REMOVED|2016] Currently kept for
#'   backward compatibility. It does not offer any speed gain. (In fact, it
#'   currently slows things down when set to anything different from
#'   \code{Inf}).
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to
#'   \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param p Penalty factor for Jaro-Winkler distance. The valid range for 
#'   \code{p} is \code{0 <= p <= 0.25}. If \code{p=0} (default), the
#'   Jaro-distance is returned. Applies only to \code{method='jw'}.
#' @param nthread Maximum number of threads to use. By default, a sensible
#'   number of threads is chosen, see \code{\link{stringdist-parallelization}}.
#'  
#' @seealso \code{\link{stringsim}}, \code{\link{qgrams}}, \code{\link{amatch}}
#'
#' @return For \code{stringdist},  a vector with string distances of size
#'   \code{max(length(a),length(b))}.
#'  
#' For \code{stringdistmatrix}: if both \code{a} and \code{b} are passed, a
#' \code{length(a)xlength(b)} \code{matrix}. If a single argument \code{a} is
#' given an object of class \code{\link[stats]{dist}} is returned.
#'  
#' Distances are nonnegative if they can be computed, \code{NA} if any of the
#' two argument strings is \code{NA} and \code{Inf} when \code{maxDist} is
#' exceeded or, in case of the hamming distance, when the two compared strings
#' have different length.
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
    warning("Argument 'maxDist' is deprecated for function 'stringdist'. This argument will be removed in the future.")   
  if (is.list(a)|is.list(b))
    warning(listwarning("stringdist","seqdist"))
            
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

  if (method == 'jw') weight <- weight[c(2,1,3)]
  do_dist(a=b, b=a
    , method=method
    , weight=weight
    , maxDist=maxDist
    , q=q
    , p=p
    , useBytes=useBytes
    , nthread=nthread)
}


#' @param useNames Use input vectors as row and column names?
#' @param ncores [DEPRECATED AND WILL BE REMOVED|2016]. Use \code{nthread} in
#'   stead. This argument is ignored.
#' @param cluster [DEPRECATED AND WILL BE REMOVED|2016].  A custom cluster,
#'   created with \code{\link[parallel]{makeCluster}}.
#'
#'
#' @rdname stringdist
#' @export
#' @rdname stringdist
stringdistmatrix <- function(a, b
  , method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard","jw","soundex")
  , useBytes = FALSE
  , weight=c(d=1,i=1,s=1,t=1),  maxDist=Inf, q=1, p=0
  , useNames=c('none','strings','names'), ncores=1, cluster=NULL
  , nthread = getOption("sd_num_thread")
){
  if (maxDist < Inf)
    warning("Argument 'maxDist' is deprecated for function 'stringdistmatrix'. This argument will be removed in the future.") 
  if (ncores > 1 ){
    warning("Argument 'ncores' is deprecated as stringdist now uses multithreading by default. This argument is currently ignored and will be removed in the future.")
    ncores <- 1
  }
  if ( !is.null(cluster) ){
    message("Argument 'cluster' is deprecaterd as stringdust now uses multithreading by default. The argument is currently ignored and will be removed in the future")
  }
  if (is.list(a)|| (!missing(b) && is.list(b)) ){
   warning(listwarning("stringdistmatrix","seqdistmatrix"))
  }

  # for backward compatability with stringdist <= 0.9.0
  if (identical(useNames, FALSE)) useNames <- "none"
  if (identical(useNames, TRUE)) useNames <- "strings"
  useNames <- match.arg(useNames)
  
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
  
  # if b is missing, generate a 'dist' object.  
  if (missing(b)){ 
    if (useNames == "names"){
      a <- setNames(as.character(a),names(a))
    } else {
      a <- as.character(a)
    }
    return( lower_tri(a
        , method=method
        , useBytes=useBytes
        , weight=weight
        , useNames=useNames
        , nthread=nthread)
    )
  }

  if (useNames == "names"){
    rowns <- names(a)
    colns <- names(b)
  }
  
  # NOTE: this strips off names
  a <- as.character(a)
  b <- as.character(b)

  if (useNames=="strings"){
    rowns <- a
    colns <- b
  } 
  
  
  if (!useBytes){
    a <- enc2utf8(a)
    b <- enc2utf8(b)
  }
 
  if (length(a) == 0 || length(b) == 0){ 
   return(matrix(numeric(0)))
  }

  x <- vapply(b, do_dist, USE.NAMES=FALSE, FUN.VALUE=numeric(length(a))
          , a, method,weight,maxDist, q, p,useBytes, nthread)

  if (useNames %in% c("strings","names") ){  
    structure(matrix(x,nrow=length(a),ncol=length(b), dimnames=list(rowns,colns)))
  } else {
    matrix(x,nrow=length(a),ncol=length(b)) 
  }
}


char2int <- function(x){
  # For some OS's enc2utf8 had unexpected behavior for NA's,
  # see https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15201.
  # This is fixed for R >= 2.15.3.
  # i <- !is.na(x)
  # x[i] <- enc2utf8(x[i])
  lapply(enc2utf8(x),utf8ToInt)
}

#  enum-type in stringdist.h
METHODS <- c(
    osa     = 0L
  , lv      = 1L
  , dl      = 2L
  , hamming = 3L
  , lcs     = 4L
  , qgram   = 5L
  , cosine  = 6L
  , jaccard = 7L
  , jw      = 8L
  , soundex = 9L
)


do_dist <- function(a, b, method, weight, maxDist=Inf, q, p, useBytes=FALSE, nthread=1L){
  
  if (method=='soundex' && !all(printable_ascii(a) & printable_ascii(b)) ){
    warning("Non-printable ascii or non-ascii characters in soundex. Results may be unreliable. See ?printable_ascii.")
  }
  method <- METHODS[method]
  if ( is.na(method) ){
    stop(sprintf("method '%s' is not defined",method))
  }

  d <- .Call("R_stringdist", a, b, method
    , as.double(weight), as.double(p), as.integer(q)
    , as.integer(useBytes), as.integer(nthread)
  )

  if (maxDist < Inf ){
    d[!is.na(d) & d > maxDist] <- Inf
  }
  d
}

# more efficient function that returns a square distance matrix as a 'stats::dist' object.
lower_tri <- function(a
   , method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard","jw","soundex")
   , useBytes = FALSE
   , weight=c(d=1,i=1,s=1,t=1), q=1, p=0
   , useNames=FALSE
   , nthread = getOption("sd_num_thread")
){
  methnr <- METHODS[method]
  if (is.na(method)){
    stop(sprintf("method '%s' is not defined",method))
  }
  
  x <- .Call("R_lower_tri", a, methnr
             , as.double(weight), as.double(p), as.integer(q)
             , as.integer(useBytes), as.integer(nthread))
  
  attributes(x) <- list(class='dist'
    , Size     = length(a)
    , Diag   = FALSE
    , Upper  = FALSE
    , method = method)
  if (useNames == "strings") attr(x,"Labels") <- as.character(a)
  if (useNames == "names" ) attr(x,"Labels") <- names(a)
  
  x
}




