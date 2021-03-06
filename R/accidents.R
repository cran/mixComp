#' Accidents Dataset
#'
#' Number of accidents incurred by 414 machinists over a period of three months
#' from Karlis and Xekalaki (1999).
#'
#' @name accidents
#' @docType data
#' @usage data(accidents)
#' @format A data frame with 414 observations on 1 variable. Replicates are generated to reflect
#' the number of accidents n incurred by michinists over a tree-month period (n = 0, 2, ..., 8).
#' As there are 296 machinists that had no accidents, 0 appears 296 times in the data,
#' as there are 74 machinists that had one accident, 1 appears 74 times in the data, etc.
#' @keywords datasets
#' @source Karlis, D., Xekalaki, E. (1999) On Testing for the Number of Components in
#' a Mixed Poisson Model. Annals of the Institute of Statistical Mathematics 51, 149-162.
#' @examples
#' data(accidents)
#'
#' \donttest{
#' # convert the data to vector:
#' accidents.obs <- unlist(accidents)
#'
#' # generate MLE function:
#' MLE.pois <- function(dat) mean(dat)
#'
#' # generate function needed for estimating the j^th moment of the
#' # mixing distribution via Hankel.method "explicit"
#'
#' explicit.pois <- function(dat, j){
#'   mat <- matrix(dat, nrow = length(dat), ncol = j) -
#'          matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
#'   return(mean(apply(mat, 1, prod)))
#' }
#'
#' # construct a 'datMix' object:
#' accidents.dM <- datMix(accidents.obs, dist = "pois", discrete = TRUE,
#'                       Hankel.method = "explicit",
#'                       Hankel.function = explicit.pois,
#'                       theta.bound.list = list(lambda = c(0, Inf)),
#'                       MLE.function = MLE.pois)
#'
#' # define the penalty:
#' pen <- function(j, n) j * log(n)
#'
#' ## complexity estimation:
#' set.seed(0)
#' res <- paramHankel(accidents.dM, j.max = 5, B = 1000, ql = 0.025, qu = 0.975)
#'
#' # plot the results:
#' plot(res, breaks = 8, ylim = c(0, 0.8))
#' }
NULL
