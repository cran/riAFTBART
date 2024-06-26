#' Perform Variable Selection using Three Threshold-based Procedures
#'
#' Performs variable selection with ri-AFTBART using the three thresholding methods introduced in Bleich et al. (2013).
#'
#' @param M.burnin A numeric value indicating the number of MCMC iterations to be treated as burn in.
#' @param M.keep A numeric value indicating the number of MCMC posterior draws after burn in.
#' @param M.thin A numeric value indicating the thinning parameter.
#' @param status A vector of event indicators: status = 1 indicates that the event was observed while status = 0 indicates the observation was right-censored.
#' @param y.train A vector of follow-up times.
#' @param x.train A dataframe or matrix, including all the covariates but not treatments for training data, with rows corresponding to observations and columns to variables.
#' @param trt.train A numeric vector representing the treatment groups for the training data.
#' @param x.test A dataframe or matrix, including all the covariates but not treatments for testing data, with  rows corresponding to observations and columns to variables.
#' @param trt.test A numeric vector representing the treatment groups for the testing data.
#' @param cluster.id A vector of integers representing the clustering id. The cluster id should be an integer and start from 1.
#' @param verbose A logical indicating whether to show the progress bar. The default is FALSE.
#' @param n_permuate Number of permutations of the event time together with the censoring indicator to generate the null permutation distribution.
#' @param alpha Cut-off level for the thresholds.
#' @param seed An optional integer seed for reproducibility of the permutations. Default is NULL.
#' 
#' @return A list with the following elements:
#' \item{var_local_selected:}{A character vector including all the variables selected using Local procedure.}
#' \item{var_max_selected:}{A character vector including all the variables selected using Global Max procedure.}
#' \item{var_global_se_selected:}{A character vector including all the variables selected using Global SE procedure.}
#' \item{vip_perm:}{The permutation distribution for the variable inclusion proportions generated by permuting the event time together with the censoring indicator.}
#' \item{vip_obs:}{The variable inclusion proportions for the actual data.}
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(20181223)
#' n = 2
#' k = 50
#' N = n*k
#' cluster.id = rep(1:n, each=k)
#' tau.error = 0.8
#' b = rnorm(n, 0, tau.error)
#' alpha = 2
#' beta1 = 1
#' beta2 = -1
#' beta3 = -2
#' sig.error = 0.5
#' censoring.rate = 0.02
#' x1 = rnorm(N,0.5,1)
#' x2 = rnorm(N,1.5,0.5)
#' error = rnorm(N,0,sig.error)
#' logtime = alpha + beta1*x1 + beta2*x2 + b[cluster.id] + error
#' y = exp(logtime)
#' C = rexp(N, rate=censoring.rate)
#' Y = pmin(y,C)
#' status = as.numeric(y<=C)
#' trt.train = sample(c(1,2,3), N, prob = c(0.4,0.3,0.2), replace = TRUE)
#' trt.test = sample(c(1,2,3), N, prob = c(0.3,0.4,0.2), replace = TRUE)
#' res <- var_select(M.burnin = 10, M.keep = 10, M.thin = 1, status = status,
#'                       y.train = Y, trt.train = trt.train, trt.test = trt.test,
#'                       x.train = cbind(x1,x2),
#'                       x.test = cbind(x1,x2),
#'                       cluster.id = cluster.id,
#'                       n_permuate = 4,alpha = 0.1,seed = 20181223)
#'                       }
var_select <- function(M.burnin, M.keep, M.thin = 1, status, y.train, x.train, trt.train, x.test, trt.test, cluster.id, verbose = FALSE, n_permuate, alpha = 0.1,seed = NULL){
  # Get the true vip from oberserved data
  riAFTBART_fit_obs <- riAFTBART_fit(M.burnin = M.burnin,
                M.keep = M.keep,
                M.thin = M.thin,
                status = status,
                y.train = y.train,
                x.train = x.train,
                trt.train = trt.train,
                x.test = x.test,
                trt.test = trt.test,
                cluster.id = cluster.id,
                verbose = verbose)
  vip_obs <- apply(riAFTBART_fit_obs$vip,2,mean) # Get the mean vip across all the posterior samples
  n_status <- length(status) # Get the length of the status
  vip_perm_matrix <- matrix(NA, nrow = n_permuate, ncol = length(vip_obs))
  for (i in 1:n_permuate){
    if (!is.null(seed)) set.seed(seed + i)
    index <- sample(1:n_status) # Get the permutation index
    riAFTBART_fit_perm <- riAFTBART_fit(M.burnin = M.burnin,
                                       M.keep = M.keep,
                                       M.thin = M.thin,
                                       status = status[index], # Permutate event time together with the censoring indicators
                                       y.train = y.train[index], # Permutate event time together with the censoring indicators
                                       x.train = x.train,
                                       trt.train = trt.train,
                                       x.test = x.test,
                                       trt.test = trt.test,
                                       cluster.id = cluster.id,
                                       verbose = verbose)
    vip_perm_matrix[i,] <- apply(riAFTBART_fit_perm$vip,2,mean)
  }
  # use local cutoff
  # alpha <- 0.1
  cutoff_local <- apply(vip_perm_matrix, 2, stats::quantile, probs = 1 - alpha)
  var_select_local = names(vip_obs[vip_obs > cutoff_local])
  # use global max cutoff
  cutoff_max = stats::quantile(apply(vip_perm_matrix, 1 ,max), 1 - alpha)
  var_select_max = names(vip_obs[vip_obs >= cutoff_max])
  # use global se cutoff
  perm_se = apply(vip_perm_matrix, 2, stats::sd)
  perm_mean = apply(vip_perm_matrix, 2, mean)
  cover_constant = bisectK(tol = .01 , coverage = 1 - alpha, permute_mat = vip_perm_matrix, x_left = 1, x_right = 20, countLimit = 100, perm_mean = perm_mean, perm_se = perm_se) # Use the method in Justin Bleich (2016) to get the constant
  var_select_global_se = names(vip_obs[which(vip_obs >= perm_mean + cover_constant * perm_se)])
  result <- list(var_local_selected = var_select_local,
                 var_max_selected = var_select_max,
                 var_global_se_selected = var_select_global_se,
                 vip_perm = vip_perm_matrix,
                 vip_obs = vip_obs)
  return(result)
}
