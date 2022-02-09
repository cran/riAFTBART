#' Fit a random effect accelerated failure time BART model
#'
#' This function is called by riAFTBART function and implement the riAFT-BART algorithm
#' @param M.burnin A numeric value indicating the number of MCMC iterations to be treated as burn in.
#' @param M.keep A numeric value indicating the number of MCMC posterior draws after burn in.
#' @param M.thin A numeric value indicating the thinning parameter.
#' @param status A vector of event indicators: status = 1 indicates that the event was observed while status = 0 indicates the observation was right-censored.
#' @param y.train A vector of follow-up times.
#' @param x.train A dataframe or matrix, including all the covariates but not treatments for training data, with rows corresponding to observations and columns to variables.
#' @param trt.train A numeric vector representing the treatment groups for the training data.
#' @param x.test A dataframe or matrix, including all the covariates but not treatments for testing data, with  rows corresponding to observations and columns to variables.
#' @param SA A logical indicating whether to conduct sensitivity analysis. The default is FALSE.
#' @param prior_c_function_used Prior confounding functions used for SA, which is inherited from the sa function. The default is NULL.
#' @param gps Generalized propensity score, which is inherited from the sa function. The default is NULL.
#' @param trt.test A numeric vector representing the treatment groups for the testing data.
#' @param cluster.id A vector of integers representing the clustering id.
#' @param verbose A logical indicating whether to show the progress bar. The default is FALSE
#' @return A list with the following elements:
#' \item{b:}{A matrix including samples from the posterior of the random effects.}
#' \item{tree:}{A matrix with M.keep rows and nrow(x.train) columns represnting the predicted log survival time for x.train.}
#' \item{tree.pred:}{A matrix with M.keep rows and nrow(x.test) columns represnting the predicted log survival time for x.test.}
#' \item{tau:}{A vector representing the posterior samples of tau, the standard deviation of the random effects.}
#' \item{alpha:}{A vector representing the posterior samples of alpha, the parameter expansion term.}
#' \item{sigma:}{A vector representing the posterior samples of sigma, the residual/error standard deviation.}
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
#'}

riAFTBART_fit <- function(M.burnin, M.keep, M.thin = 1, status, y.train, x.train, trt.train, x.test, trt.test, cluster.id, verbose = FALSE, SA = FALSE,prior_c_function_used = NULL, gps = NULL)
{
  # initial values
  N = length(y.train)  # total sample size
  n = length(unique(cluster.id)) # number of clusters
  null_aft_model <- survival::survreg(survival::Surv(y.train, status) ~ 1, dist="lognormal")
  null_intercept <- null_aft_model$coefficients
  aft_model_scale <- null_aft_model$scale
  aft_model <- survival::survreg(survival::Surv(y.train, status) ~ ., dist="lognormal",data = data.frame(y.train, x.train))
  b.init = tapply(stats::resid(aft_model), cluster.id, mean)
  sigma.init <- stats::sd(tapply(stats::resid(aft_model), cluster.id, mean))
  x.train <- cbind(x.train, trt.train)
  colnames(x.train) <- NULL
  x.test <- cbind(x.test, trt.test)
  colnames(x.test) <- NULL
  tau.init = 1
  alpha.init = 1
  tree.init = stats::rnorm(N)
  logT.init = log(y.train)-null_intercept
  # logT.init = log(y.train)
  c.id <- which(status==0)        # censored.id
  n.censored <- length(c.id)
  logT.censored <- NULL
  for(c in 1:n.censored) logT.censored[c] <- msm::rtnorm(1, mean=logT.init[c.id[c]], sd=sigma.init, lower=logT.init[c.id[c]], upper=Inf)
  logT.init[status==0] <- logT.censored

  # Define and initialize the list of things to keep track of in the "current state" of the chain
  cur <- list(tree=tree.init, b=b.init, sigma=sigma.init, tau=tau.init, alpha = alpha.init, logT=logT.init)

  # Define matrix of store MCMC results
  # T <- data$T
  # X <- data$X
  P <- ncol(x.train)        # number of covariates
  N.pred <- nrow(x.test) # sample size of test dataset

  chain.keep <- list()
  chain.keep$b <- array(NA, dim = c(n, M.keep))       # random effects
  chain.keep$tree <- array(NA, dim = c(N, M.keep))    # tree structure
  chain.keep$tree.pred <- array(NA, dim = c(N.pred, M.keep))    # predicted tree structure
  chain.keep$tau <- array(NA, dim = c(1,M.keep))      # standard deviation of random intercept
  chain.keep$alpha <- array(NA, dim = c(1,M.keep))      # parameter expansion term
  chain.keep$sigma <- array(NA, dim = c(1,M.keep))    # standard deviation of residual error
  #chain.keep$logT <- array(NA, dim = c(N,M.keep))     # log-transformation of event times, including censored times

  rownames(chain.keep$tau) <- 'tau'
  rownames(chain.keep$alpha) <- 'alpha'
  rownames(chain.keep$sigma) <- 'sigma'
  rownames(chain.keep$b) <- paste('b', 1:n, sep="")
  rownames(chain.keep$tree) <- paste('tree', 1:N, sep="")
  rownames(chain.keep$tree.pred) <- paste('tree.pred', 1:N.pred, sep="")
  #rownames(chain.keep$logT) <- paste('Unit', 1:N, sep="")
  if (verbose == TRUE){
  cat("Running the MCMC algorithm. It may be long, keep cool :)", "\n")
  cat("\n")
}
  # Burn-in phase: do not keep these results
  if (M.burnin > 0) {
    if (verbose == TRUE){
    cat("Burn-in phase", "\n")
    }
    for (i in 1:M.burnin) {
      cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,aft_model_scale,SA = SA,trt.train,prior_c_function_used,gps)
      if (verbose == TRUE){
        progressBar(i, M.burnin)
      }

    } # end of i for-loop
  } # end of if-loop

  # Converged phase: keep these results after thinning
  if (verbose == TRUE){
  cat("\n")
  cat("Converged phase", "\n")
}
  for (m in 1:M.keep) {
    if (M.thin > 1) {	# Skip the "thinned" pieces of the chain
      for (j in 1:(M.thin-1)) {
        cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,aft_model_scale,SA = SA,trt.train,prior_c_function_used,gps)
      }
    }
    cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,aft_model_scale,SA = SA,trt.train,prior_c_function_used,gps)
    if (verbose == TRUE){
      progressBar(m, M.keep)
    }
    chain.keep$b[,m] <- cur$b
    chain.keep$tree[,m] <- cur$tree
    chain.keep$tau[m] <- cur$tau
    chain.keep$alpha[m] <- cur$alpha
    chain.keep$sigma[m] <- cur$sigma
    #chain.keep$logT[,m] <- cur$logT
    chain.keep$tree.pred[,m] <- cur$tree.pred

  } # end of m for-loop
  chain.keep$tau <- c(chain.keep$tau)
  chain.keep$sigma <- c(chain.keep$sigma)
  chain.keep$alpha <- c(chain.keep$alpha)
  chain.keep$tree <- chain.keep$tree + null_intercept
  chain.keep$tree.pred <- chain.keep$tree.pred + null_intercept
  class(chain.keep) <- "riAFTBART_estimate"
  return(chain.keep)

}
