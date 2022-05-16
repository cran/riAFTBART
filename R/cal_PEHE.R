#' Calculate the PEHE
#'
#' This function calculates the PEHE based on the survival probability from a fitted ri-AFTBART model.
#'
#' @param object An object from cal_survprob() function.
#' @param metric A character string representing the metric to be calculated for PEHE. Only \code{"survival"} and \code{"rmst"} is allowed.
#' @param time A numeric value representing the time point used to calculate PEHE.
#' @param LP A numeric vector corresponding to the true linear predictors for each treatment from the simulated data.
#' @param lambda A numeric value representing the true follow up time for from the simulated data.
#' @param eta A numeric value to induce proportional/non-proportional hazards assumption from the simulated data.
#'
#' @return A list with the following three components:
#' \item{true:}{A numeric vector representing the true survival or rmst for each individual.}
#' \item{predicted:}{A numeric vector representing the predicted survival or rmst for each individual.}
#' \item{pehe:}{A numeric vector representing the calculated pehe.}
#' @export
#'
#' @examples
#' \donttest{
#' library(riAFTBART)
#' lp_w_all <-
#'   c(".4*x1 + .1*x2  - .1*x4 + .1*x5",    #' w = 1
#'     ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5")  #' w = 2
#' nlp_w_all <-
#'   c("-.5*x1*x4  - .1*x2*x5", #' w = 1
#'     "-.3*x1*x4 + .2*x2*x5")#' w = 2
#' lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
#' nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#' X_all <- c(
#'   "rnorm(10, 0, 0.5)",#' x1
#'   "rbeta(10, 2, .4)",   #' x2
#'   "runif(10, 0, 0.5)",#' x3
#'   "rweibull(10,1,2)",  #' x4
#'   "rbinom(10, 1, .4)"#' x5
#' )
#' set.seed(111111)
#' data <- dat_sim(
#'   nK = 2,
#'   K = 5,
#'   n_trt = 3,
#'   X = X_all,
#'   eta = 2,
#'   lp_y = lp_y_all,
#'   nlp_y  = nlp_y_all,
#'   align = FALSE,
#'   lp_w = lp_w_all,
#'   nlp_w = nlp_w_all,
#'   lambda = c(1000,2000,3000),
#'   delta = c(0.5,0.5),
#'   psi = 1,
#'   sigma_w = 1,
#'   sigma_y = 2,
#'   censor_rate = 0.1
#' )
#' data$LP_true[,1]
#' data$lambda
#' data$eta
#' res <- riAFTBART_fit(M.burnin = 10, M.keep = 10, M.thin = 1, status = data$delta,
#'                       y.train = data$Tobs, trt.train = data$w, trt.test = 1,
#'                       x.train = data$covariates,
#'                       x.test = data$covariates,
#'                       cluster.id = data$cluster)
#' res_cal_surv_prob <- cal_surv_prob(object = res,
#' time.points = 1:max(data$Tobs),
#' test.only = TRUE,
#' cluster.id = data$cluster)
#'
#' res_cal_PEHE_survival <- cal_PEHE(object = res_cal_surv_prob,
#'                          metric = "survival", time = 40,
#'                          LP = data$LP_true[,1], lambda = data$lambda[1],
#'                          eta = data$eta)
#'
#' res_cal_PEHE_rmst <- cal_PEHE(object = res_cal_surv_prob,
#'                                   metric = "rmst",
#'                                   time = 40,
#'                                   LP = data$LP_true[,1],
#'                                   lambda = data$lambda[1],
#'                                   eta = data$eta)
#'                                   }

cal_PEHE <- function(object, metric, time, LP, lambda, eta){
  fun_true <-
    function(x) {
      exp(-(1/lambda * exp(LP[i]) * x) ^ eta)
    }
  if (metric == "survival"){
    true <- NULL
    predicted <- NULL
    pehe <- NULL
    for (i in 1:length(LP)){
      true[i] <- fun_true(time)
      predicted[i] <- object$Surv.test[i,][object$time.points == time]
      pehe[i] <- true[i]-predicted[i]
    }
    return(list(true = true, predicted = predicted, pehe = pehe))
  } else if (metric == "rmst"){
    true <- NULL
    predicted <- NULL
    pehe <- NULL
    for (i in 1:length(LP)){
      integrate_res <- stats::integrate(fun_true, lower = 0, upper = time)
      true[i] <- integrate_res$value
      predicted[i] <- RISCA::rmst(object$time.points, object$Surv.test[i,], max.time = time, type = "l")
      pehe[i] <- true[i]-predicted[i]
    }
    return(list(true = true, predicted = predicted, pehe = pehe))
  }
}
