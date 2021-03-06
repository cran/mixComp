#' Children Dataset
#'
#' Number of children of 4075 widows entitled to support from a certain pension fund from Thisted (1988).
#'
#' @name children
#' @docType data
#' @usage data(children)
#' @format A data frame with 4075 observations on 1 variable. Replicates are generated
#' to reflect the number of children n that widows entitled to support have  (n = 0, 1, ..., 6).
#' As there are 3062 widows that have no children, 0 appears 3062 times in the data,
#' as there are 587 widows that have one child, 1 appears 587 times in the data, etc.
#' @keywords datasets
#' @source Thisted, R. A. (1988). Elements of statistical computing:
#' Numerical computation (Vol. 1). CRC Press.
#' @examples
#' data(children)
#'
#' # convert the data to vector:
#' children.obs <- unlist(children)
#' 
#' # explicit function giving the estimate for the j^th moment of the
#' # mixing distribution, needed for Hankel.method "explicit"
#' explicit.pois <- function(dat, j){
#'   mat <- matrix(dat, nrow = length(dat), ncol = j) -
#'          matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
#'   return(mean(apply(mat, 1, prod)))
#' }
#' 
#' # define the MLE function:
#' MLE.pois <- function(dat) mean(dat)
#' 
#' # construct a 'datMix' object:
#' children.dM <- datMix(children.obs, dist = "pois", discrete = TRUE,
#'                       Hankel.method = "explicit",
#'                       Hankel.function = explicit.pois,
#'                       theta.bound.list = list(lambda = c(0, Inf)),
#'                       MLE.function = MLE.pois)
#'                       
#' # define the penalty:
#' pen <- function(j, n) j * log(n)
#' 
#' # complexity estimation:
#' \donttest{
#' set.seed(0)
#' det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE,
#'                                B = 1000, pen.function = pen)
#' plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset")
#' }
NULL
