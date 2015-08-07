#' @title 
#' String metrics in \pkg{stringdist}
#' 
#' @description
#' This page gives an overview of the string dissimilarity measures offered by
#' \pkg{stringdist}.
#' 
#' @section String Metrics:
#' String metrics are ways of quantifying the dissimilarity between two finite
#' sequences, usually text strings. Over the years, many such measures have been
#' developed. Some are based on a mathematical understanding of the set of all 
#' strings that can be composed from a finite alphabet, others are based on more
#' heuristic principles, such as how a text string sounds when pronounced by a 
#' native English speaker.
#'
#' The terms 'string metrics' and 'string distance' are used more or less
#' interchangibly in literature. From a mathematical point of view, string
#' metrics often do not obey the demands that are usually required from a
#' distance function. For example, it is not true for all string metrics that a
#' distance of 0 means that two strings are the same (e.g. in the \eqn{q}-gram
#' distance). Nevertheless, string metrics are very useful in practice and have
#' many applications.
#'
#' The metric you need to choose for an application strongly depends on both the
#' nature of the string (what does the string represent?) and the cause of
#' dissimilarities between the strings you are measuring. For example, if you
#' are comparing human-typed names that may contain typo's, the Jaro-Winkler
#' distance may be of use. If you are comparing names that were written down
#' after hearing them, a phonetic distance may be a better choice.
#'
#' Currently, the following distance metrics are supported by \pkg{stringdist}.
#' \tabular{ll}{
#'    \bold{Method name} \tab \bold{Description}\cr
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
#'
#' @section A short description of string metrics supported by \pkg{stringdist}:
#'
#' See \href{http://journal.r-project.org/archive/2014-1/loo.pdf}{Van der Loo
#' (2014)} for an extensive description and references. The review papers of
#' Navarro (2001) and Boytsov (2011) provide excellent technical overviews of
#' respectively online and offline string matching algorithms.
#' 
#' The \bold{Hamming distance} (\code{method='hamming'}) counts the number of 
#' character substitutions that turns \code{b} into \code{a}. If \code{a} 
#' and \code{b} have different number of characters the distance is \code{Inf}. 
#'
#' The \bold{Levenshtein distance} (\code{method='lv'}) counts the number of 
#' deletions, insertions and substitutions necessary to turn \code{b} into 
#' \code{a}. This method is equivalent to \code{R}'s native \code{\link[utils]{adist}} 
#' function. 
#'
#' The \bold{Optimal String Alignment distance} (\code{method='osa'}) is like the Levenshtein 
#' distance but also allows transposition of adjacent characters. Here, each 
#' substring  may be edited only once. (For example, a character cannot be transposed twice
#' to move it forward in the string). 
#'
#' The \bold{full Damerau-Levensthein distance} (\code{method='dl'}) is like the optimal 
#' string alignment distance except that it allows for multiple edits on substrings. 
#'
#' The \bold{longest common substring} (method='lcs') is defined as the longest string that can be 
#' obtained by pairing characters from \code{a} and \code{b} while keeping the order 
#' of characters intact. The \bold{lcs-distance} is defined as the number of unpaired characters. 
#' The distance is equivalent to the edit distance allowing only deletions and insertions, 
#' each with weight one. 
#'
#' A \bold{\eqn{q}-gram} (method='qgram') is a subsequence of \eqn{q} \emph{consecutive} 
#' characters of a string. If \eqn{x} (\eqn{y}) is the vector of counts
#' of \eqn{q}-gram occurrences in \code{a} (\code{b}), the \bold{\eqn{q}-gram distance} 
#' is given by the sum over the absolute differences \eqn{|x_i-y_i|}.
#' The computation is aborted when \code{q} is is larger than the length of 
#' any of the strings. In that case \code{Inf}  is returned.
#'
#' The \bold{cosine distance} (method='cosine') is computed as \eqn{1-x\cdot y/(\|x\|\|y\|)}, where \eqn{x} and 
#' \eqn{y} were defined above.
#' 
#' Let \eqn{X} be the set of unique \eqn{q}-grams in \code{a} and \eqn{Y} the set of unique 
#' \eqn{q}-grams in \code{b}. The \bold{Jaccard distance} (\code{method='jaccard'}) is given by \eqn{1-|X\cap Y|/|X\cup Y|}.
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
#' correction term to the Jaro-distance. It is defined as \eqn{d - l\cdot p\cdot d}, where
#' \eqn{d} is the Jaro-distance. Here,  \eqn{l} is obtained by counting, from
#' the start of the input strings, after how many characters the first
#' character mismatch between the two strings occurs, with a maximum of four. The
#' factor \eqn{p} is a penalty factor, which in the work of Winkler is often
#' chosen \eqn{0.1}.
#'
#' For the \bold{soundex} distance (method='soundex'), strings are translated to a soundex code 
#' (see \code{\link{phonetic}} for a specification). The
#' distance between strings is 0 when they have the same soundex code,
#' otherwise 1. Note that soundex recoding is only meaningful for characters
#' in the ranges a-z and A-Z. A warning is emitted when non-printable or non-ascii
#' characters are encountered. Also see \code{\link{printable_ascii}}.
#'
#'
#'
#' @references
#'
#' \itemize{
#' \item{
#'   MPJ van der Loo (2014) \emph{The stringdist package for approximate string matching}. The R Journal \bold{6}(1) 111-122.
#' }
#' \item{
#'  L. Boytsov (2011). \emph{Indexing methods for approximate dictionary searching: comparative analyses}. ACM Journal of experimental
#'  algorithmics \bold{16} 1-88.
#' }
#' \item{
#'  G. Navarro (2001). \emph{A guided tour to approximate string matching}. ACM Computing Surveys \bold{33} 31-88.
#' }
#' }
#' @seealso
#' \itemize{ 
#'  \item{Functions applying string metrics to text: \code{\link{stringdist}},
#'    \code{\link{stringdistmatrix}}, \code{\link{amatch}}}
#'  \item{Functions applying string metrics to integer sequences:
#'    \code{\link{seq_dist}}, \code{\link{seq_distmatrix}}, \code{\link{seq_amatch}} }
#'  \item{Encoding issues: \code{\link{stringdist-encoding}}  }
#' }
#' 
#' @name stringdist-metrics
#' @rdname stringdist-metrics
{}


