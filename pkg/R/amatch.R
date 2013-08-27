#' Approximate string matching
#'
#' Approximate string matching equivalents of \code{R}'s native \code{\link[base]{match}} and \code{\%in\%}.
#'
#' @section Note on \code{NA} handling:
#' \code{R}'s native \code{\link[base]{match}} function matches \code{NA} with \code{NA}. This
#' may feel inconsistent with \code{R}'s usual \code{NA} handling, since for example \code{NA==NA} yields
#' \code{NA} rather than \code{TRUE}. In most cases, one may reason about the behaviour
#' under \code{NA} along the lines of ``if one of the arguments is \code{NA}, the result shall be \code{NA}'',
#' simply because not all information necessary to execute the function is available. One uses special 
#' functions such as \code{is.na}, \code{is.null} \emph{etc.} to handle special values. 
#'
#' The \code{amatch} function mimics the behaviour of \code{\link[base]{match}} by default: \code{NA} is 
#' matched with \code{NA} and with nothing else. Note that this is inconsistent with the behaviour of \code{\link{stringdist}}
#' since \code{stringdist} yields \code{NA} when at least one of the arguments is \code{NA}. The same inconsistency exists
#' between \code{\link[base]{match}} and \code{\link[stats]{dist}}. However, in \code{amatch} this behaviour
#' can be controlled by setting \code{matchNA=FALSE}. In that case, if any of the arguments in \code{x} 
#' is \code{NA}, the \code{nomatch} value is returned, regardless of whether \code{NA} is present in \code{table}.
#'
#'
#' @param x vector: elements to be approximately matched
#' @param table vector: lookup table for matching
#' @param nomatch the value to be returned when no match is found. This is coerced to integer. \code{nomatch=0} 
#'  can be a useful option.
#' @param matchNA Should \code{NA}'s be matched? Default behaviour mimics the
#'   behaviour of base \code{\link[base]{match}}, meaning that \code{NA} matches
#'   \code{NA} (which is inconsistent with \code{dist} or \code{stringdist}).
#' @param method Matching algorithm to use. See \code{\link{stringdist}}.
#' @param useBytes Perform byte-wise comparison. \code{useBytes=TRUE} is faster but may yield different
#' 	results depending on character encoding. For \code{ASCII} it is identical. See also \code{\link{strindist}}.
#' @param weight parameters for matching algorithm See \code{\link{stringdist}}.
#' @param maxDist Elements in \code{x} will not be matched with elements of
#'  \code{table} if their distance is larger than \code{maxDist}. 
#'   
#' @param q q-gram size, see \code{\link{stringdist}}.
#' @param p Winklers penalty parameter for Jaro-Winkler distance, see \code{\link{stringdist}}.
#' @return \code{amatch} returns the position of the closest match of \code{x} in \code{table}. 
#'  When multiple matches with the same smallest distance metric exist, the first one is returned.
#'  \code{ain} returns a \code{logical} vector of length \code{length(x)} indicating wether 
#'  an element of \code{x} approximately matches an element in \code{table}.
#'
#' @example ../examples/amatch.R
#' @export
amatch <- function(x, table, nomatch=NA_integer_, matchNA=TRUE, 
  method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard", "jw"), 
  useBytes = FALSE,
  weight=c(d=1,i=1,s=1,t=1), 
  maxDist=0.1, q=1, p=0){


  x <- as.character(x)
  table <- as.character(table)

  method <- match.arg(method)
  if (!useBytes){
    x <- char2int(x)
    table <- char2int(table)
  }
  stopifnot(
      all(is.finite(weight))
      , all(weight > 0)
      , all(weight <=1)
      , q >= 0
      , p <= 0.25
      , p >= 0
      , matchNA %in% c(TRUE,FALSE)
      , maxDist > 0
      , is.logical(useBytes)
  )
  if (maxDist==Inf) maxDist <- 0L;
  switch(method,
    osa     = .Call('R_match_osa'       , x, table, as.integer(nomatch), as.integer(matchNA), as.double(weight), as.double(maxDist)),
    lv      = .Call('R_match_lv'        , x, table, as.integer(nomatch), as.integer(matchNA), as.double(weight), as.double(maxDist)),
    dl      = .Call('R_match_dl'        , x, table, as.integer(nomatch), as.integer(matchNA), as.double(weight), as.double(maxDist)),
    hamming = .Call('R_match_hm'        , x, table, as.integer(nomatch), as.integer(matchNA), as.integer(maxDist)),
    lcs     = .Call('R_match_lcs'        , x, table, as.integer(nomatch), as.integer(matchNA), as.integer(maxDist)),
    qgram   = .Call('R_match_qgram_tree', x, table, as.integer(nomatch), as.integer(matchNA), as.integer(q), as.double(maxDist), 0L),
    cosine  = .Call('R_match_qgram_tree', x, table, as.integer(nomatch), as.integer(matchNA), as.integer(q), as.double(maxDist), 1L),
    jaccard = .Call('R_match_qgram_tree', x, table, as.integer(nomatch), as.integer(matchNA), as.integer(q), as.double(maxDist), 2L),
    jw      = .Call('R_match_jw'        , x, table, as.integer(nomatch), as.integer(matchNA), as.double(p), as.double(maxDist))
  )
}

#' @param ... parameters to pass to \code{amatch} (except \code{nomatch})
#' @rdname amatch
#' @export 
ain <- function(x,table,...){
  amatch(x, table, nomatch=0, ...) > 0
}

