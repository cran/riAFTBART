bisectK <- function(tol, coverage, permute_mat, x_left, x_right, countLimit, perm_mean, perm_se){
  count = 0
  guess = mean(c(x_left, x_right))
  while ((x_right - x_left) / 2 >= tol & count < countLimit){
    empirical_coverage = mean(sapply(1 : nrow(permute_mat),
                                     function(s){all(permute_mat[s,] - perm_mean <= guess * perm_se)}))
    if (empirical_coverage - coverage == 0){
      break
    } else if (empirical_coverage - coverage < 0){
      x_left = guess
    } else {
      x_right = guess
    }
    guess = mean(c(x_left, x_right))
    count = count + 1
  }
  guess
}
