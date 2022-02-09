#' Plot the propensity score by treatment
#'
#' This function estimates the propensity score for each treatment group and then plot the propensity score by each treatment to check covariate overlap.
#'
#' @param trt A numeric vector representing the treatment groups.
#' @param X A dataframe or matrix, including all the covariates but not treatments, with  rows corresponding to observations and columns to variables.
#' @param method A character indicating how to estimate the propensity score. The default is "Multinomial", which uses multinomial regression to estimate the propensity score.
#'
#' @return A plot
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
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
#' plot_gps(trt = trt.train, X = cbind(x1, x2))
plot_gps <- function(trt, X, method = "Multinomial"){
  group <- NULL
  ps <- NULL
  es.max.ATE <- NULL
  if (method == "Multinomial"){
    multinom_result <- nnet::multinom(trt ~ ., data =   as.data.frame(cbind(trt, X)))
    pred_ps <- stats::fitted(multinom_result)
    p_1 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 1) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 1|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p_2 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 2) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 2|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p_3 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 3) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 3|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    (p <- cowplot::plot_grid(p_1, p_2, p_3,ncol = 3,align = "h"))
    return(p)
  } else if (method == "GBM"){
    X <- as.data.frame(X)
    temp<- noquote(names(X))
    strFormula  = sprintf("trt~%s", paste(temp, sep = "",collapse="+"))
    psmod <- twang::mnps(stats::as.formula(strFormula),
                         data=as.data.frame(cbind(trt, X)) %>% dplyr::mutate(trt = as.factor(trt)), estimand = "ATE")

    p_1 <- psmod$psList$`1`$ps %>%
        cbind(trt) %>%
        dplyr::mutate(trt = as.factor(trt)) %>%
        ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
        ggplot2::geom_boxplot()+
        ggplot2::labs(x = "", y = "", title = "P(A = 1|X, V)")+
        ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
        ggplot2::theme_bw()+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p_2 <- psmod$psList$`2`$ps %>%
      cbind(trt) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 2|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p_3 <- psmod$psList$`3`$ps %>%
      cbind(trt) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 3|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    (p <- cowplot::plot_grid(p_1, p_2, p_3,ncol = 3,align = "h"))
    return(p)
  }
}
