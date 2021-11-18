#' @describeIn forward_probabilities Matrix of backwards log-probabilities
backward_logprobabilities <- function(x, Gamma, p){
  #p_mat <- outer(1:ncol(Gamma), x, FUN=p)
  p_mat <- matrix(0, nrow=nrow(Gamma), ncol=length(x))
  for(i in 1:nrow(Gamma)){
    for(j in 1:length(x)){
      p_mat[i, j] <- p(i, x[j])
    }
  }
  return(backward_ll_cpp(log(Gamma), log(p_mat)))
}
