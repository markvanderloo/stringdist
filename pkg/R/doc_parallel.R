#' @title
#' Multithreading and parallelization in \pkg{stringdist}
#'
#' 
#' @description This page describes how \pkg{stringdist} uses parallel processing.
#'
#' @section Multithreading and parallelization in \pkg{stringdist}: The core 
#'   functions of \pkg{stringdist} are implemented in C. On systems where 
#'   \code{openMP} is available, \pkg{stringdist} will automatically take 
#'   advantage of multiple cores. The
#'   \href{https://cran.r-project.org/doc/manuals/r-release/R-exts.html#OpenMP-support}{section
#'   on OpenMP} of the
#'   \href{https://cran.r-project.org/doc/manuals/r-release/R-exts.html}{Writing
#'   R Extensions} manual discusses on what systems OpenMP is available (at the time of writing more or
#'   less, anywhere except on OSX).
#'
#' By default, the number of threads to use is taken from \code{options('sd_num_thread')}.
#' When the package is loaded, the value for this option is determined as follows:
#' \itemize{
#'   \item{If the environment variable \code{OMP_NUM_THREADS} is set, this value is taken.}
#'   \item{Otherwise, the number of available cores is determined with \code{parallel::detectCores()}
#'   If this fails, the number of threads is set to 1 (with a message). If the nr of detected
#'   cores exceeds three, the number of used cores is set to \eqn{n-1}.}
#'   \item{If available, the environment variable \code{OMP_THREAD_LIMIT} is
#'   determined and The number of threads is set to the lesser of
#'   \code{OMP_THREAD_LIMIT} and the number of detected cores.}
#' }
#'
#' The latter step makes sure that on machines with \eqn{n>3} cores, \eqn{n-1} 
#' cores are used. Some benchmarking showed that using all cores is often slower
#' in such cases. This is probably because at least one of the threads will be
#' shared with the operating system.
#'
#' Functions that use multithreading have an option named \code{nthread} that
#' controls the maximum number of threads to use. If you need to do large
#' calculations, it is probably a good idea to benchmark the performance on your
#' machine(s) as a function of \code{'nthread'}, for example using the 
#' \href{https://cran.r-project.org/package=microbenchmark}{microbenchmark}
#' package of Mersmann.
#'
#'
#'
#'
#' @seealso
#' \itemize{
#'  \item{Functions running multithreaded: \code{\link{stringdist}}, \code{\link{stringdistmatrix}}, \code{\link{amatch}}, \code{\link{ain}} }
#' }
#'
#' @name stringdist-parallelization
#' @rdname stringdist-parallelization
{}
