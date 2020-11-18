## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mixComp)

## -----------------------------------------------------------------------------
set.seed(0)
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))

## ---- fig.width = 5, fig.height = 4-------------------------------------------
plot(normLocMix, main = "Three component normal mixture")
plot(poisMix, main = "Three component poisson mixture")

## -----------------------------------------------------------------------------
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
# plot the histograms of the random samples:
plot(normLocRMix, main = "Three component normal mixture")
plot(poisRMix, main = "Three component poisson mixture")

## -----------------------------------------------------------------------------
obs <- faithful$waiting
norm.dist <- "norm"

## -----------------------------------------------------------------------------
norm.bound.list <- vector(mode = "list", length = 2)
names(norm.bound.list) <- c("mean", "sd")
norm.bound.list$mean <- c(-Inf, Inf)
norm.bound.list$sd <- c(0, Inf)

## -----------------------------------------------------------------------------
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)

## -----------------------------------------------------------------------------
method <- "translation"
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

## -----------------------------------------------------------------------------
faithful.dM <- datMix(obs, dist = norm.dist, 
                      theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)

## -----------------------------------------------------------------------------
explicit.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}

## -----------------------------------------------------------------------------
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}

## -----------------------------------------------------------------------------
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

## -----------------------------------------------------------------------------
set.seed(1)
geomMix <- Mix("geom", w = c(0.1, 0.6, 0.3), prob = c(0.8, 0.2, 0.4))
geomRMix <- rMix(1000, obj = geomMix)
geom.dM <- RtoDat(geomRMix, Hankel.method = "explicit", 
                  Hankel.function = explicit.geom)

## -----------------------------------------------------------------------------
pen <- function(j, n){
  (j * log(n)) / (sqrt(n))
}
set.seed(1)
geomdets_sca_pen <- nonparamHankel(geom.dM, j.max = 5, scaled = TRUE, 
                                   B = 1000, pen.function = pen)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
print(geomdets_sca_pen)
plot(geomdets_sca_pen, main = "Three component geometric mixture")

## -----------------------------------------------------------------------------
MLE.geom <- function(dat){
  1/(mean(dat)+1)
}
geom.dM <- RtoDat(geomRMix, Hankel.method = "explicit", 
                  Hankel.function = explicit.geom, 
                  theta.bound.list = list(prob = c(0, 1)), 
                  MLE.function = MLE.geom)
set.seed(1)
res <- paramHankel(geom.dM, j.max = 5, B = 1000, ql = 0.025, qu = 0.975)

## -----------------------------------------------------------------------------
normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list)
set.seed(0)
res <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000, 
                      threshold = "AIC")

## ---- fig.width = 5, fig.height = 4-------------------------------------------
print(res)
plot(res)

## ---- results='hide', message=FALSE, warning=FALSE----------------------------
poisList <- vector(mode = "list", length = 1)
names(poisList) <- "lambda"
poisList$lambda <- c(0, Inf)
MLE.pois <- function(dat) mean(dat)
pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, 
                  MLE.function = MLE.pois)
set.seed(1)
res <- hellinger.boot.disc(pois.dM, B = 50, ql = 0.025, qu = 0.975)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
print(res)
plot(res)

## ---- fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE----
set.seed(1)
res <- mix.lrt(faithful.dM, B = 50, quantile = 0.95)
print(res)
plot(res)

## ---- fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE----
dnorm0.5 <- function(x, mean){
  dnorm(x, mean = mean,  sd = 0.5)
}
rnorm0.5 <- function(n, mean){
  rnorm(n, mean = mean,  sd = 0.5)
}
## create objects `Mix` and `rMix`:
set.seed(1)
norm0.5Mix <- Mix("norm0.5", w = c(0.3, 0.4, 0.3), mean = c(10, 11, 13))
norm0.5RMix <- rMix(1000, obj = norm0.5Mix)
## plot the results:
plot(norm0.5Mix)
plot(norm0.5RMix)

## -----------------------------------------------------------------------------
norm0.5.list <- vector(mode = "list", length = 1)
names(norm0.5.list) <- c("mean")
norm0.5.list$mean <- c(-Inf, Inf)

MLE.norm0.5 <- function(dat) mean(dat)

norm0.5.dM <- RtoDat(norm0.5RMix, theta.bound.list = norm0.5.list,
                     MLE.function = MLE.norm0.5)

## ---- results='hide', message=FALSE, warning=FALSE----------------------------
set.seed(1)
res <- mix.lrt(norm0.5.dM, B = 50, quantile = 0.95)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
print(res)
plot(res)

