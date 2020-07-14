#' Stringdist-based fuzzy text search
#'
#' \code{afind} slides a window of fixed width over a string \code{x} and
#' computes the distance between the current window and the sought-after
#' \code{pattern}. The location, content, and distance corresponding to the
#' window with the best match is returned.
#'
#'
#' @param x \code{[character]} strings to search in
#' @param pattern \code{[character]} strings to find (not a regular expression).
#' @param window \code{[integer]} width of moving window
#' @param value \code{[logical]} toggle return matrix with matched strings.
#' @inheritParams amatch
#'
#' @details
#' Matching is case-sensitive.  Both \code{x} and \code{pattern} are converted
#' to \code{UTF-8} prior to search, unless \code{useBytes=TRUE}, in which case
#' the distances are measured bytewise.
#'
#' Code is parallelized over the \code{x} variable: each value of \code{x}
#' is scanned for every element in \code{pattern} using a separate thread (when \code{nthread}
#' is larger then 1).
#'
#' The current implementation is naive, in the sense that for each string
#' \code{s} in \code{x}, \code{nchar(s) - window + 1} separate distances are
#' computed. At the moment no attempt is made to speed up the calculation by
#' using that consecutive windows overlap.
#'
#'
#' @return
#' A \code{list} of three matrices, each of with \code{length(x)} rows and \code{length(pattern)}
#' columns. In each matrix, element \eqn{(i,j)} corresponds to \code{x[i]} and \code{pattern[j]}.
#' \itemize{
#' \item{\code{location}. \code{[integer]}, location of the start of best matching window.
#'       When \code{useBytes=FALSE}, this corresponds to the location of a \code{UTF} code point
#'       in \code{x}, possibly after conversion from its original encoding.}
#' \item{\code{distance}. \code{[character]}, the string distance between pattern and
#'       the best matching window.}
#' \item{\code{match}. \code{[character]}, the first, best matching window.}
#' 
#' }
#' 
#' @family matching
#'
#' @examples
#' texts = c("When I grow up, I want to be"
#'        , "one of the harversters of the sea"
#'        , "I think before my days are gone"
#'        , "I want to be a fisherman")
#' patterns = c("fish", "gone","to be")
#'
#' afind(texts, patterns, method="cosine", q=3)
#'
#'
#' @export
afind <- function(x, pattern, window=nchar(enc2utf8(pattern))
  , value=TRUE
  , method = c("osa","lv","dl","hamming","lcs", "qgram","cosine","jaccard","jw","soundex")
  , useBytes = FALSE
  , weight=c(d=1,i=1,s=1,t=1) 
  , q  = 1
  , p  = 0
  , bt = 0
  , nthread = getOption("sd_num_thread")
  ){
  
  stopifnot(
    all(is.finite(weight))
    , all(weight > 0)
    , all(weight <=1)
    , window > 0
    , q >= 0
    , p <= 0.25
    , p >= 0
    , is.logical(useBytes)
    , ifelse(method %in% c('osa','dl'), length(weight) >= 4, TRUE)
    , ifelse(method %in% c('lv','jw') , length(weight) >= 3, TRUE)
    , nthread > 0
  )
  x <- as.character(x)
  pattern <- as.character(pattern)
  if ( !useBytes ){
    x <- enc2utf8(x)
    pattern <- enc2utf8(pattern)
  }

  if (length(x) == 0) return(numeric(0))

  method <- match.arg(method)
  if (method == 'jw') weight <- weight[c(2,1,3)]

  method <- METHODS[method]
  if ( is.na(method) ){
    stop(sprintf("method '%s' is not defined",method))
  }
  
  L <- .Call("R_afind"
    , x
    , pattern
    , as.integer(window)
    , method
    , as.double(weight)
    , as.double(p)
    , as.double(bt)
    , as.integer(q)
    , as.integer(useBytes)
    , as.integer(nthread)
    , PACKAGE="stringdist")
  
  names(L) <- c("location", "distance")

  if (isTRUE(value)){
    matches = sapply(seq_along(pattern), function(i){ 
      substr(x, L[[1]][,i], L[[1]][,i] + window[i]-1)
    })
    L$match <- matrix(matches, nrow=length(x))
  }

  L
}



