#' Plot the Trace plots for the paramters from a fitted riAFT-BART model
#'
#' This function create the trace plots for the parameters from a fitted riAFT-BART model.
#'
#' @param x A fitted object of from riAFTBART_estimate function.
#' @param focus A character specifying which parameter to plot.
#' @param id A numeric vector indicating the subject or cluster index to plot
#' @param ... further arguments passed to or from other methods.
#'
#' @return A plot
#' @export
#'
#' @examples
#' \donttest{
#' library(riAFTBART)
#' set.seed(20181223)
#' n = 50      # number of clusters
#' k = 50      # cluster size
#' N = n*k     # total sample size
#' cluster.id = rep(1:n, each=k)

#' tau.error = 0.8
#' b = stats::rnorm(n, 0, tau.error)

#' alpha = 2
#' beta1 = 1
#' beta2 = -1
#' sig.error = 0.5
#' censoring.rate = 0.02

#' x1 = stats::rnorm(N,0.5,1)
#' x2 = stats::rnorm(N,1.5,0.5)
#' trt.train = sample(c(1,2,3), N, prob = c(0.4,0.3,0.2), replace = TRUE)
#' trt.test = sample(c(1,2,3), N, prob = c(0.3,0.4,0.2), replace = TRUE)
#' error = stats::rnorm(N,0,sig.error)
#' logtime = alpha + beta1*x1 + beta2*x2 + b[cluster.id] + error
#' y = exp(logtime)
#' C = rexp(N, rate=censoring.rate) # censoring times
#' Y = pmin(y,C)
#' status = as.numeric(y<=C)
#' res <- riAFTBART_fit(M.burnin = 50, M.keep = 50, M.thin = 1, status = status,
#'                       y.train = Y, trt.train = trt.train, trt.test = trt.test,
#'                       x.train = cbind(x1,x2),
#'                       x.test = cbind(x1,x2),
#'                       cluster.id = cluster.id)
#' plot(x = res, focus = "sigma")
#' }
plot.riAFTBART_estimate <- function(x, focus = "sigma", id = NULL,...){
  if (focus == "sigma"){
    plot(x$sigma, type="l", main="sigma", ylab="value", col="gray")
  } else if (focus == "tau"){
    plot(x$tau, type="l", main="tau", ylab="value", col="gray")
  } else if (focus == "alpha"){
    plot(x$alpha, type="l", main="alpha", ylab="value", col="gray")
  } else if (focus == "b"){
    plot(x$b[id,], type="l", main="Random intercept", ylab="value", col="gray")
  } else if (focus == "tree"){
    plot(x$tree[id,], type="l", main="Tree (train)", ylab="value", col="gray")
  } else if (focus == "tree.pred"){
    plot(x$tree[id,], type="l", main="Tree (test)", ylab="value", col="gray")
  }

}
