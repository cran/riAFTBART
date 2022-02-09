#' Flexible Monte Carlo sensitivity analysis for unmeasured confounding
#'
#' The function \code{sa} implements the flexible sensitivity analysis approach for unmeasured confounding with multiple treatments from multilevel survival data.
#' @param M.burnin A numeric value indicating the number of MCMC iterations to be treated as burn in.
#' @param M.keep A numeric value indicating the number of MCMC posterior draws after burn in.
#' @param M.thin A numeric value indicating the thinning parameter.
#' @param status A vector of event indicators: status = 1 indicates that the event was observed while status = 0 indicates the observation was right-censored.
#' @param y.train A vector of follow-up times.
#' @param x.train A dataframe or matrix, including all the covariates but not treatments for training data, with rows corresponding to observations and columns to variables.
#' @param trt.train A numeric vector representing the treatment groups for the training data.
#' @param formula A \code{\link[stats]{formula}} object for the analysis. The default is to use all terms specified in \code{x.train}.
#' @param x.test A dataframe, including all the covariates but not treatments for testing data, with rows corresponding to observations and columns to variables.
#' @param trt.test A numeric vector representing the treatment groups for the testing data.
#' @param cluster.id A vector of integers representing the clustering id.
#' @param verbose A logical indicating whether to show the progress bar. The default is FALSE
#' @param prior_c_function 1) A vector of characters indicating the prior distributions for the confounding functions. Each character contains the random number generation code from the standard probability \code{\link[stats:Distributions]{distributions}} in the \code{\link[stats:stats-package]{stats}} package. 2) A vector of characters including the grid specifications for the confounding functions. It should be used when users want to formulate the  confounding  functions as scalar values. 3) A matrix indicating the point mass prior for the confounding functions
#' @param Q1 A numeric value indicating the number of draws of the GPS from the posterior predictive distribution
#' @param Q2 A numeric value indicating the number of draws from the prior distributions of the confounding functions
#' @param nCores A numeric value indicating number of cores to use for parallel computing.
#' @param ... Other parameters that can be passed to BART functions
#'
#' @return A list of causal estimands including risk difference in terms of log time between different treatment groups.
#'
#' @export
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' set.seed(20181223)
#' n = 50      # number of clusters
#' k = 50      # cluster size
#' N = n*k     # total sample size
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
#' trt.train = sample(c(1,2,3), N, prob = c(0.4,0.3,0.2), replace = TRUE)
#' trt.test = sample(c(1,2,3), N, prob = c(0.3,0.4,0.2), replace = TRUE)
#' error = rnorm(N,0,sig.error)
#' logtime = alpha + beta1*x1 + beta2*x2 + b[cluster.id] + error
#' y = exp(logtime)
#' C = rexp(N, rate=censoring.rate) # censoring times
#' Y = pmin(y,C)
#' status = as.numeric(y<=C)
#' res_sa <- sa(M.burnin = 50, M.keep = 50, M.thin = 1, status = status,
#'              y.train = Y,trt.train = trt.train,trt.test = trt.test,
#'              x.train = cbind(x1,x2),
#'              x.test = cbind(x1,x2),
#'              cluster.id = cluster.id, verbose = F,prior_c_function = c(
#'                "runif(-0.6, 0)",# c(1,2)
#'                "runif(0, 0.6)",# c(2,1)
#'                "runif(-0.6, 0)", # c(2,3)
#'                "seq(-0.6, 0, by = 0.3)", # c(1,3)
#'                "seq(0, 0.6, by = 0.3)", # c(3,1)
#'               "runif(0, 0.6)" # c(3,2)
#'             ),Q1 = 1, nCores = 1)
#'  }
sa <- function(M.burnin, M.keep, M.thin = 1, status, y.train, x.train, trt.train, x.test, trt.test, cluster.id, verbose = FALSE, formula = NULL, prior_c_function, Q1, Q2 = NULL, nCores = 1,... ){
  if (sum(c(length(trt.train) == length(y.train), length(trt.train) == nrow(x.train), length(y.train) == nrow(x.train))) != 3) stop(paste0("The length of y.train, the length of trt.train and the nrow for x.train should be equal. Please double check the input."), call. = FALSE)
  if (!is.null(formula)){
    x.train <- as.data.frame(stats::model.matrix(object = formula, cbind(y.train,x.train)))
    x.train <- x.train[,!(names(x.train) == "(Intercept)")]
  }
  if (any(stringr::str_detect(prior_c_function, "seq")) == FALSE && is.numeric(prior_c_function) == FALSE){
    prior_c_function_all <- matrix(NA, ncol = length(prior_c_function), nrow = Q2)
    for (i in 1:length(prior_c_function)){
      str_locate_parenthesis <- stringr::str_locate(prior_c_function[i],"\\(")
      # set.seed(seed)
      prior_c_function_all[,i] <- eval(parse(text = paste0(paste0(stringr::str_sub(prior_c_function[i], 1, str_locate_parenthesis[1]), Q2, ",",stringr::str_sub(prior_c_function[i], str_locate_parenthesis[1]+1)))))
    }
    prior_c_function_used <- prior_c_function_all
  }
  if (any(stringr::str_detect(prior_c_function, "seq")) == TRUE){
    c_index_with_grid <- which(stringr::str_detect(prior_c_function, "seq"))
    c_index_without_grid <- which(!stringr::str_detect(prior_c_function, "seq"))
    n_c_with_grid <- length(prior_c_function[stringr::str_detect(prior_c_function, "seq")])
    grid_length <- length(eval(parse(text = prior_c_function[stringr::str_detect(prior_c_function, "seq")])))
    Q2 <- grid_length^n_c_with_grid
    c_with_grid <- prior_c_function[c_index_with_grid]
    c_without_grid <- prior_c_function[c_index_without_grid]
    c_without_grid_all <- matrix(NA, ncol = length(c_without_grid), nrow = Q2)
    for (i in 1:length(c_without_grid)){
      str_locate_parenthesis <- stringr::str_locate(c_without_grid[i],"\\(")
      # set.seed(seed)
      c_without_grid_all[,i] <- eval(parse(text = paste0(paste0(stringr::str_sub(c_without_grid[i], 1, str_locate_parenthesis[1]), Q2, ",",stringr::str_sub(c_without_grid[i], str_locate_parenthesis[1]+1)))))
    }
    colnames(c_without_grid_all) <- c_index_without_grid
    c_with_grid_1 <- NULL
    for (i in 1:length(c_index_with_grid)){
      assign(paste0("c_with_grid_",i), eval(parse(text = c_with_grid[i])))
    }
    c_with_grid_all <- c_with_grid_1
    for (i in 1:(length(c_index_with_grid)-1)){
      c_with_grid_all <- tidyr::expand_grid(c_with_grid_all, eval(parse(text = paste0("c_with_grid_",(i+1)))))
    }
    colnames(c_with_grid_all) <- c_index_with_grid
    c_functions_grid_final <- cbind(as.data.frame(c_without_grid_all), c_with_grid_all)
    names(c_functions_grid_final) <- paste0("c", names(c_functions_grid_final) )
    c_functions_grid_final <- c_functions_grid_final %>%
      dplyr::select(paste0("c",1:length(prior_c_function)))
    prior_c_function_used <- c_functions_grid_final
  }
  if (all(stringr::str_detect(prior_c_function, "runif"))) {
    prior_c_function_used <- matrix(NA, ncol = length(prior_c_function), nrow = Q2)
    for (i in 1:length(prior_c_function)){
      str_locate_parenthesis <- stringr::str_locate(prior_c_function[i],"\\(")
      # set.seed(seed)
      prior_c_function_used[,i] <- eval(parse(text = paste0(paste0(stringr::str_sub(prior_c_function[i], 1, str_locate_parenthesis[1]), Q2, ",",stringr::str_sub(prior_c_function[i], str_locate_parenthesis[1]+1)))))
    }
  }

  if (is.numeric(prior_c_function) == TRUE){
    prior_c_function_used <- t(apply(prior_c_function, 2, mean))
  }
  # change the type of y.train and trt.train as the input parameter of bart function

  x.train_A_model = cbind(as.data.frame(x.train), cluster.id) %>% dplyr::mutate_if(is.character, as.factor)
  y.train = as.numeric(y.train)
  log.y.train <- log(y.train)
  trt.train = as.integer(trt.train)
  n_trt <- length(unique(trt.train))
  prior_c_function_used = as.matrix(prior_c_function_used)
  n_alpha = nrow(prior_c_function_used)
  # set.seed(seed)
  # fit the treatment assigment model, to use gap-sampling, we over sample n * 10 samples, and select a sample per 10 turns
  A_model <- BART::mbart2(x.train = x.train_A_model, as.integer(as.factor(trt.train)), x.test = x.train_A_model, ndpost = Q1 * 10, mc.cores = nCores)
  # assign the estimated assignment probability to each sample, the size is (n, #treatment, sample_size)
  gps = array(A_model$prob.test[seq(1, nrow(A_model$prob.test), 10),], dim = c(Q1, length(unique(trt.train)), length(trt.train)))
  train_x = cbind(x.train, trt.train)
  n_trt <- length(unique(trt.train))
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  iterations <- n_alpha
    out <-
      foreach::foreach(
        i = 1:n_alpha,
        .combine = function(x, y) {
          list(result_riAFTBART = rbind(x[["result_riAFTBART"]], y[["result_riAFTBART"]]))
        }
      ) %dopar% {

        # cat('Starting ', i, 'th job.\n', sep = '')
        result_riAFTBART <- NULL
        for (j in 1:Q1) {
          # fit the bart model to estimate causal effect
          bart_mod = riAFTBART::riAFTBART_fit(M.burnin = M.burnin, M.keep = M.keep, M.thin = 1, status = status,
                                    y.train = exp(log.y.train), trt.train = as.numeric(trt.train), trt.test = as.numeric(trt.test),
                                    x.train = x.train,
                                    x.test = x.test,
                                    cluster.id = cluster.id, verbose = F,
                                    SA = TRUE,
                                    gps = gps[j,,],
                                    prior_c_function_used = prior_c_function_used[i,])
          result_riAFTBART_once <- rowMeans(bart_mod$tree.pred)
          names(result_riAFTBART_once) <- NULL
          result_riAFTBART <- c(result_riAFTBART, result_riAFTBART_once)
        }
        result_riAFTBART = list(result_riAFTBART = result_riAFTBART)
      }
    parallel::stopCluster(cl)
    result_riAFTBART <- out$result_riAFTBART
    if (any(stringr::str_detect(prior_c_function, "seq")) == TRUE){
      counter <- 1
      result_riAFTBART <- list(result_riAFTBART = result_riAFTBART, c_functions = prior_c_function_used, grid_index = c_index_with_grid)

      return(result_riAFTBART)
    }
    if (any(stringr::str_detect(prior_c_function, "seq")) == FALSE){

      return(result_riAFTBART)
    }
}

