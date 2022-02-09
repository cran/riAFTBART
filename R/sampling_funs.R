# update tree structure, BART(X) and residual standard error
sample.tree <- function(status, y.train, x.train, x.test, cluster.id, old,aft_model_scale, trt.train,SA = SA,prior_c_function_used,gps) {
  logT <- old$logT
  b <- old$b

  # fit BART using Z as the outcome
  Z <- logT - b[cluster.id]
  if (SA){ # Sensitivity analysis to correct for Z
    train_x = cbind(x.train, trt.train)
    n_trt <- length(unique(trt.train))
    Z = ifelse(train_x[, "trt.train"] == sort(unique(train_x[, "trt.train"]))[1], Z - (unlist(prior_c_function_used[1]) * gps[2, ] + unlist(prior_c_function_used[4]) * gps[3, ]),
                         ifelse(train_x[, "trt.train"] == sort(unique(train_x[, "trt.train"]))[2], Z - (unlist(prior_c_function_used[2]) * gps[1, ] + unlist(prior_c_function_used[3]) * gps[3, ]),
                                Z - (unlist(prior_c_function_used[5]) * gps[1, ] + unlist(prior_c_function_used[6]) * gps[2, ])))
  }

  nskip <- 5
  ndpost <- 2
  mod <- dbarts::bart(x.train=x.train, y.train=Z, x.test=x.test, ntree=200, ndpost=ndpost, nskip=nskip, keepevery=1, verbose=FALSE,k=aft_model_scale)

  tree <- mod$yhat.train[ndpost,]
  sigma <- mod$sigma[ndpost]
  tree.pred <- mod$yhat.test[ndpost,]

  cur <- old
  cur$tree <- tree
  cur$sigma <- sigma
  cur$tree.pred <- tree.pred
  cur$SA <- SA
  return(cur)
}



# update random intercept, b
sample.b <- function(status, y.train, x.train, x.test, cluster.id, old) {

  # X <- data$X
  # cluster.id <- data$cluster.id
  n <- length(unique(cluster.id))
  cluster.size <- as.numeric(table(cluster.id))

  logT <- old$logT
  tau2 <- old$tau^2
  alpha <- old$alpha
  sigma2 <- old$sigma^2
  tree <- old$tree

  b.new <- NULL
  for (k in 1:n) {
    b.var <- (cluster.size[k]/sigma2 + 1/(tau2 * alpha))^(-1)
    yy <- logT[cluster.id==k]
    tt <- tree[cluster.id==k]
    resid <- yy - tt
    b.mu <- b.var * sum(resid) / sigma2
    b.new[k] <- stats::rnorm(1, b.mu, sqrt(b.var))
  }

  cur <- old
  cur$b <- b.new
  return(cur)
}

# update parameter expansion term
# Inverse-Gamma prior, alpha ~ IG(1,1)
# alpha <- NULL
sample.alpha <- function(status, y.train, x.train, x.test, cluster.id, old) {

  b <- old$b
  tau2 <- old$tau^2
  # cluster.id <- data$cluster.id
  n <- length(unique(cluster.id))


  rss <- t(b) %*% b
  alpha <- MCMCpack::rinvgamma(1, shape = 1, scale = rss/(2 * tau2) + 1)

  cur <- old
  cur$alpha <- sqrt(alpha)
  return(cur)
}



# update standard deviation of random intercept, tau
# Inverse-Gamma prior, tau2 ~ IG(d1,d2)
sample.tau <- function(status, y.train, x.train, x.test, cluster.id, old) {
  hyperpars = list(d1=1, d2=1)
  b <- old$b
  # cluster.id <- data$cluster.id
  n <- length(unique(cluster.id))
  alpha <- old$alpha
  d1 <- hyperpars$d1
  d2 <- hyperpars$d2

  rss <- t(b) %*% b
  tau2 <- MCMCpack::rinvgamma(1, shape = n/2 + d1, scale = rss/(2*alpha) + d2)
  # tau2 <- MCMCpack::rinvgamma(1, shape = n/2 + d1, scale = rss/2 + d2)
  cur <- old
  cur$tau <- sqrt(tau2)
  return(cur)
}



# impute censored event time
sample.time.censored <- function(status, y.train, x.train, cluster.id, old) {

  # T <- data$T                     # event time
  # Status <- data$Status           # event status
  lower.bound <- log(y.train)
  c.id <- which(status==0)        # censored.id
  n.censored <- length(c.id)
  # cluster.id <- data$cluster.id

  # X <- data$X
  N <- nrow(x.train)

  b <- old$b
  tree <- old$tree
  sigma <- old$sigma
  logT <- old$logT

  logT.mean <- tree + b[cluster.id]

  logT.censored <- NULL
  for(c in 1:n.censored) logT.censored[c] <- msm::rtnorm(1, mean=logT.mean[c.id[c]], sd=sigma, lower=lower.bound[c.id[c]], upper=Inf)
  logT[status==0] <- logT.censored

  cur <- old
  cur$logT <- logT
  return(cur)
}
