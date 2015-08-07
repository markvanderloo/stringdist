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
#' \code{\link{stringdist}} and then dividing the distance by the maximum
#' possible distance. This results in a score between 0 and 1, with 1
#' corresponding to perfect similarity and 0 to complete dissimilarity.
#' For distances where weights can be specified, the maximum distance is currently computed by
#' assuming that all weights are equal to 1.
#'
#' @example ../examples/stringsim.R
#' @export
stringsim <- function(a, b, method = c("osa", "lv", "dl", "hamming", "lcs",
  "qgram", "cosine", "jaccard", "jw", "soundex"), useBytes=FALSE, q = 1, ...) {
  # Calculate the distance 
  method <- match.arg(method)
  dist <- stringdist::stringdist(a, b, method=method, useBytes=useBytes, q=q, ...)

  nctype <- if (useBytes) "bytes" else "char"

  # Normalise the distance by dividing it by the maximum possible distance
  if (method == "hamming") {
    max_dist <- if (length(b) > length(a)) nchar(b,type=nctype) else nchar(a,type=nctype)
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "lcs") {
    max_dist <- nchar(a,type=nctype) + nchar(b,type=nctype)
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "lv") {
    max_dist <- pmax(nchar(a,type=nctype), nchar(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "osa") {
    max_dist <- pmax(nchar(a,type=nctype), nchar(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "dl") {
    max_dist <- pmax(nchar(a,type=nctype), nchar(b,type=nctype))
    max_dist[max_dist == 0] <- 1
    sim <- 1 - dist/max_dist
  } else if (method == "qgram") {
    max_dist <- (nchar(a,type=nctype) + nchar(b,type=nctype) - 2*q + 2)
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

