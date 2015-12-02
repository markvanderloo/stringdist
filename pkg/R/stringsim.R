#' Compute similarity scores between strings
#'
#' \code{stringsim} computes pairwise string similarities between elements of
#' \code{character} vectors \code{a} and \code{b}, where the vector with less
#' elements is recycled. 
#'
#' @param a R object (target); will be converted by \code{as.character}.
#' @param b R object (source); will be converted by \code{as.character}.
#' @param method Method for distance calculation. The default is \code{"osa"}, 
#'   see \code{\link{stringdist-metrics}}.
#' @param useBytes Perform byte-wise comparison, see \code{\link{stringdist-encoding}}.
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to
#'   \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param ... additional arguments are passed on to \code{\link{stringdist}}.
#'
#' @return
#' Returns a vector with similarities, which are values between 0 and 1 where
#' 1 corresponds to perfect similarity (distance 0) and 0 to complete
#' dissimilarity. \code{NA} is returned when \code{\link{stringdist}} returns 
#' \code{NA}. Distances equal to \code{Inf} are truncated to a similarity of
#' 0.
#'
#' @details
#' The similarity is calculated by first calculating the distance using
#' \code{\link{stringdist}}, dividing the distance by the maximum
#' possible distance, and substracting the result from 1. 
#' This results in a score between 0 and 1, with 1
#' corresponding to complete similarity and 0 to complete dissimilarity.
#' Note that complete similarity only means equality for distances satisfying
#' the identity property. This is not the case e.g. for q-gram based distances
#' (for example if q=1, anagrams are completely similar).
#' For distances where weights can be specified, the maximum distance 
#' is currently computed by assuming that all weights are equal to 1.
#'
#' @example ../examples/stringsim.R
#' @export
stringsim <- function(a, b, method = c("osa", "lv", "dl", "hamming", "lcs",
  "qgram", "cosine", "jaccard", "jw", "soundex"), useBytes=FALSE, q = 1, ...) {
  # Calculate the distance 
  method <- match.arg(method)
  dist <- stringdist::stringdist(a, b, method=method, useBytes=useBytes, q=q, ...)

  nctype <- if (useBytes) "bytes" else "char"
  normalize_dist(dist, a, b, method=method, nctype=nctype, q=q)
}


#' Compute similarity scores between sequences of integers
#' 
#' @param a \code{list} of \code{integer} vectors (target)
#' @param b \code{list} of \code{integer} vectors (source). Optional for
#'   \code{seq_distmatrix}.
#' @param method Method for distance calculation. The default is \code{"osa"}, 
#'   see \code{\link{stringdist-metrics}}.
#' @param q  Size of the \eqn{q}-gram; must be nonnegative. Only applies to
#'   \code{method='qgram'}, \code{'jaccard'} or \code{'cosine'}.
#' @param ... additional arguments are passed on to \code{\link{seq_dist}}.
#'
#' @return 
#' A \code{numeric} vector of length \code{max(length(a),length(b))}. If one of the
#' entries in \code{a} or \code{b} is \code{NA_integer_}, all comparisons with that
#' element result in \code{NA}. Missings occurring within the sequences are treated
#' as an ordinary number (the representation of \code{NA_integer_}).
#'   
#' @example ../examples/seq_sim.R
#' @seealso \code{\link{seq_dist}}, \code{\link{seq_amatch}}
#' @export
seq_sim <- function(a, b, method = c("osa", "lv", "dl", "hamming", "lcs",
   "qgram", "cosine", "jaccard", "jw"),  q = 1, ...) {
    
  method <- match.arg(method)
  dist <- stringdist::seq_dist(a, b, method=method, q=q, ...)
  normalize_dist(dist,a,b,method=method,q=q)
}


#### HELPER FUNCTIONS ---------------------------------------------------------

# get lengths of sequences (internal function)
lengths <- function(x,...){
  UseMethod("lengths")
}

lengths.character <- function(x, type="char",...){
  nchar(x,type=type)
}

lengths.list <- function(x,...){
  .Call("R_lengths",x)
}

normalize_dist <- function(dist, a, b, method, nctype="char",q=1L){

  # Normalise the distance by dividing it by the maximum possible distance
  if (method == "hamming") {
    max_dist <- if (length(b) > length(a)) lengths(b,type=nctype) else lengths(a,type=nctype)
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "lcs") {
    max_dist <- lengths(a,type=nctype) + lengths(b,type=nctype)
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "lv") {
    max_dist <- pmax(lengths(a,type=nctype), lengths(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "osa") {
    max_dist <- pmax(lengths(a,type=nctype), lengths(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "dl") {
    max_dist <- pmax(lengths(a,type=nctype), lengths(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "qgram") {
    max_dist <- (lengths(a,type=nctype) + lengths(b,type=nctype) - 2*q + 2)
    max_dist[max_dist < 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "cosine") {
    sim <- 1 - dist
  } else if (method == "jaccard") {
    sim <- 1 - dist
  } else if (method == "jw") {
    sim <- 1 - dist
  } else if (method == "soundex") {
    sim <- 1 - dist
  }
  # all metrics can have distances == Inf; for similariy score set these to 0
  sim[sim < 0] <- 0
  sim
}



