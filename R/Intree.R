#' Interpreting Tree Ensembles with inTrees
#'
#' The inTrees (interpretable trees) framework that extracts, measures, prunes and selects rules from a tree ensemble. All the codes we use are from the inTrees github repository to act as a work around method since package inTrees was removed from the CRAN repository.
#'
#' @param X A matrix indicating the predictor variables.
#' @param Y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param ntree Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param typeDecay An integer of 1 or 2. 1 representing relative error and 2 representing error. The default is set to 2.
#' @param digits An integer indicating the digits for rounding in Intrees.
#' @param n_rule An integer indicating the minimum number of rules to consider in Intrees.
#'
#' @return A matrix including a set of relevant and non-redundant rules, and their metrics
#' @export
#'
#' @examples
#'
#' X <- within(iris,rm("Species")); Y <- iris[,"Species"]
#' intree_result <- intree(X, Y, ntree=100, digits = 3, n_rule = 2000)
#'
#'
intree <- function(X, Y, ntree, typeDecay = 2,digits, n_rule){
  rf <- randomForest::randomForest(X, Y, ntree = ntree)
  tree_list <- RF2List(rf)
  rule_exec <- extractRules(tree_list,X,digits=digits)
  rule_exec <- unique(rule_exec) # remove same rules. NOTE: for variable interaction analysis, you should NOT perform this step
  ix <- sample(1:length(rule_exec),min(n_rule,length(rule_exec))) #randomly select n_rule
  rule_exec <- rule_exec[ix,,drop=FALSE]
  rule_metric <- getRuleMetric(rule_exec,X,Y)
  rule_metric <- pruneRule(rule_metric,X,Y,typeDecay = typeDecay)
  rule_select <- selectRuleRRF(rule_metric, X, Y)
}
