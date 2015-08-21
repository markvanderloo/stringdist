#' Compute distance metrics between integer sequences
#'
#' \code{seq_dist} computes pairwise string distances between elements of 
#' \code{a} and \code{b}, where the argument with less elements is recycled. 
#' \code{seq_distmatrix} computes the distance matrix with rows according to
#' \code{a} and columns according to \code{b}.
#'
#'
#' @section Notes:
#' Input vectors are converted with \code{as.integer}. This causes truncation for numeric
#' vectors (e.g. \code{pi} will be treated as \code{3L}).
#'
#' @param a (\code{list} of) \code{integer} or \code{numeric} vector(s). Will be converted with \code{as.integer}  (target)
#' @param b (\code{list} of) \code{integer} or \code{numeric} vector(s). Will be converted with \code{as.integer} (source). 
#'    Optional for \code{seq_distmatrix}.
#' @param method Distance metric. See \code{\link{stringdist-metrics}}
#' @param weight For \code{method='osa'} or \code{'dl'}, the penalty for
#'   deletion, insertion, substitution and transposition, in that order. When
#'   \code{method='lv'}, the penalty for transposition is ignored. When
#'   \code{method='jw'}, the weights associated with characters of \code{a},
#'   characters from \code{b} and the transposition weight, in that order. 
#'   Weights must be positive and not exceed 1. \code{weight} is ignored
#'   completely when \code{method='hamming'}, \code{'qgram'}, \code{'cosine'},
#'   \code{'Jaccard'}, or \code{'lcs'}
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to
#'   \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param p Penalty factor for Jaro-Winkler distance. The valid range for 
#'   \code{p} is \code{0 <= p <= 0.25}. If \code{p=0} (default), the
#'   Jaro-distance is returned. Applies only to \code{method='jw'}.
#' @param nthread Maximum number of threads to use. By default, a sensible
#'   number of threads is chosen, see \code{\link{stringdist-parallelization}}.
#'
#' @return 
#' 
#' \code{seq_dist} returns a numeric vector with pairwise distances between \code{a}
#' and \code{b} of length \code{max(length(a),length(b)}.
#' 
#' For \code{seq_distmatrix} there are two options. If \code{b} is missing, the 
#' \code{\link[stats]{dist}} object corresponding to the \code{length(a) X
#' length(a)} distance matrix is returned. If \code{b} is specified, the
#' \code{length(a) X length(b)} distance matrix is returned.
#'    
#' If any element of \code{a} or \code{b} is \code{NA_integer_}, the distance with
#' any matched integer vector will result in \code{NA}. Missing values in the sequences
#' themselves are treated as a number and not treated specially (Also see the examples).
#'   
#' @seealso \code{\link{seq_sim}}, \code{\link{seq_amatch}}, \code{\link{seq_qgrams}} 
#'  
#' @example ../examples/seq_dist.R
#' @export
seq_dist <- function(a, b
  , method=c("osa","lv","dl","hamming","lcs", "qgram","cosine","jaccard","jw")
  , weight=c(d=1,i=1,s=1,t=1) 
  , q=1, p=0
  , nthread = getOption("sd_num_thread")
){
  a <- ensure_int_list(a)
  b <- ensure_int_list(b)
  
  stopifnot(
     all(is.finite(weight))
    , all(weight > 0)
    , all(weight <=1)
    , q >= 0
    , p <= 0.25
    , p >= 0
    , ifelse(method %in% c('osa','dl'), length(weight) >= 4, TRUE)
    , ifelse(method %in% c('lv','jw') , length(weight) >= 3, TRUE)
    , nthread > 0
  )
  
  
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
    , q=q
    , p=p
    , nthread=nthread)
}

#' @param useNames label the output matrix with \code{names(a)} and \code{names(b)}?
#' @rdname seq_dist
#' @export 
seq_distmatrix <- function(a, b
   , method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard","jw")
   , weight=c(d=1,i=1,s=1,t=1),  q=1, p=0
   , useNames=c("names","none")
   , nthread = getOption("sd_num_thread")
){
  useNames <- match.arg(useNames)
  method <- match.arg(method)
  nthread <- as.integer(nthread)
  if (method == 'jw') weight <- weight[c(2,1,3)]
  
  a <- ensure_int_list(a)
  
  # if b is missing, generate a 'dist' object.  
  if (missing(b)){ 
    return( lower_tri(a
        , method=method
        , weight=weight
        , useNames=useNames
        , nthread=nthread)
    )
  }
 
  b <- ensure_int_list(b)
  if (length(a) == 0 || length(b) == 0){ 
    return(matrix(numeric(0)))
  }
  
  if (useNames == "names"){
    rowns <- names(a)
    colns <- names(b)
  }
  
  
  #x <- vapply(b, do_dist, USE.NAMES=FALSE, FUN.VALUE=numeric(length(a))
  #      , b=a, method=method, weight=weight, q=q, p=p, nthread=nthread)
  
  x <- vapply(b
      , function(src) do_dist(list(src), b=a, method=method, weight=weight, q=q, p=p, nthread=nthread)
      , USE.NAMES=FALSE, FUN.VALUE=numeric(length(a))
    )
  
  if (useNames == "names" ){  
    structure(matrix(x,nrow=length(a),ncol=length(b), dimnames=list(rowns,colns)))
  } else {
    matrix(x,nrow=length(a),ncol=length(b)) 
  }
  
}

