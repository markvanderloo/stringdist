#'
#' @param x vector: elements to be approximately matched
#' @param table vector: lookup table for fuzzy matching
#' @param nomatch the value to be returned when no match is found. This is coerced to integer. \code{nomatch=0} 
#'  can be a useful option.
#' @param matchNA Determines wheter NA's are to be matched as in R's \code{match} function.
#' @rdname stringdist
amatch <- function(x, table, nomatch=NA_integer_, matchNA=TRUE, 
  method=c("osa","lv","dl","hamming","lcs","qgram","cosine","jaccard", "jw"), 
  weight=c(d=1,i=1,s=1,t=1), 
  maxDist=Inf, q=1, p=0){


  x <- as.character(x)
  table <- as.character(table)

  method <- match.arg(method)
  x <- char2int(x)
  table <- char2int(table)
  stopifnot(
      all(is.finite(weight)),
      all(weight > 0),
      all(weight <=1),
      p <= 0.25,
      p >= 0,
      matchNA %in% c(TRUE,FALSE)
  )
  if (maxDist==Inf) maxDist <- 0L;
  switch(method,
    osa     = .Call('R_match_osa'   , x, table,as.integer(nomatch), as.integer(matchNA), as.double(weight), as.double(maxDist))
#    lv      = .Call('R_lv'    , a, b, as.double(weight), as.double(maxDist)),
#    dl      = .Call('R_dl'    , a, b, as.double(weight), as.double(maxDist)),
#    hamming = .Call('R_hm'    , a, b, as.integer(maxDist)),
#    lcs     = .Call('R_lcs'   , a, b, as.integer(maxDist)),
#    qgram   = .Call('R_qgram_tree' , a, b, as.integer(q), 0L),
#    cosine  = .Call('R_qgram_tree' , a, b, as.integer(q), 1L),
#    jaccard = .Call('R_qgram_tree' , a, b, as.integer(q), 2L),
#    jw      = .Call('R_jaro_winkler', a, b, as.double(p))
  )
}

