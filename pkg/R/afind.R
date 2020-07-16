#' Stringdist-based fuzzy text search
#'
#' \code{afind} slides a window of fixed width over a string \code{x} and
#' computes the distance between the each window and the sought-after
#' \code{pattern}. The location, content, and distance corresponding to the
#' window with the best match is returned.
#'
#'
#' @param x  strings to search in
#' @param pattern strings to find (not a regular expression). For \code{grab},
#' \code{grabl}, and \code{extract} this must be a single string.
#' @param window  width of moving window.
#' @param value toggle return matrix with matched strings.
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
#' The current implementation of all distances except for
#' \code{"running_cosine"} is naive, in the sense that for each string \code{s}
#' in \code{x}, \code{nchar(s) - window + 1} separate distances are computed.
#' At the moment no attempt is made to speed up the calculation by using that
#' consecutive windows overlap.
#'
#' The functions \code{grab} and \code{grabl} are approximate string matching
#' functions that mimic base R's \code{\link[base]{grep}} and
#' \code{\link[base:grep]{grepl}}. They are implemented as convenience wrappers
#' of \code{find}.
#'
#' @section Running cosine distance:
#' This algorithm gains efficiency by using that two consecutive windows have
#' a large overlap in their q-gram profiles. It gives the same result as
#' the \code{"cosine"} distance, but much faster.
#'
#'
#' @return
#' For \code{afind}: a \code{list} of three matrices, each with
#' \code{length(x)} rows and \code{length(pattern)} columns. In each matrix,
#' element \eqn{(i,j)} corresponds to \code{x[i]} and \code{pattern[j]}. The 
#' names and description of each matrix is as follows.
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
#'        , "one of the harvesters of the sea"
#'        , "I think before my days are gone"
#'        , "I want to be a fisherman")
#' patterns = c("fish", "gone","to be")
#'
#' afind(texts, patterns, method="running_cosine", q=3)
#'
#' grabl(texts,"grew", maxDist=1)
#' extract(texts, "harvested", maxDist=3)
#'
#'
#' @export
afind <- function(x, pattern, window=NULL
  , value=TRUE
  , method = c("osa","lv","dl","hamming","lcs", "qgram","cosine","running_cosine","jaccard","jw","soundex")
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
    , is.null(window) || window >= 1
    , q >= 0
    , p <= 0.25
    , p >= 0
    , is.logical(useBytes) && !is.na(useBytes)
    , is.logical(value) && !is.na(value)
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

  if (is.null(window)){ 
    window = nchar(pattern, type = if (useBytes) "bytes" else "char")
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




#' @rdname afind
#' @param ... passed to \code{afind}.
#' @param maxDist Only windows with distance \code{<= maxDist} are considered a match.
#' @return 
#' For \code{grab}, an \code{integer} vector, indicating in which elements of
#' \code{x} a match was found with a distance \code{<= maxDist}. The matched
#' values when \code{value=TRUE} (equivalent to \code{\link[base]{grep}}).
#' @export
grab <- function(x, pattern, maxDist=Inf, value=FALSE, ...){
  stopifnot(is.numeric(maxDist), maxDist >= 0, length(pattern) == 1)
  L <- afind(x, pattern, value=value, ...)
  if (!value){
    which(L$distance <= maxDist)
  } else {
    L$match[L$distance <= maxDist ]
  }
}

#' @rdname afind
#' @param ... passed to \code{afind}.
#' @return 
#' For \code{grabl}, a \code{logical} vector, indicating in which elements of
#' \code{x} a match was found with a distance \code{<= maxDist}.  (equivalent
#' to \code{\link[base:grep]{grepl}}).
#' @export
grabl <- function(x, pattern, maxDist=Inf, ...){
  stopifnot(is.numeric(maxDist), maxDist >= 0, length(pattern) == 1)
  L <- afind(x, pattern, value=FALSE, ...)
  L$distance <= maxDist
}


#' @rdname afind
#'
#' @return
#' For \code{extract}, a \code{character} vector of \code{length(x)}, with
#' \code{NA} where no match was found and the first matched string if there is
#' a match. (similar to \code{stringr::str_extract}).
#' @export
extract <- function(x, pattern, maxDist = Inf, ...){
  stopifnot(is.numeric(maxDist), maxDist >= 0, length(pattern) == 1)
  L <- afind(x, pattern, value=TRUE, ...)
  out <- L$match
  out[L$distance > maxDist] <- NA_character_
  out
}


