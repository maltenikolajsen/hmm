#' @describeIn forward_probabilities Matrix of forward log-probabilities
forward_logprobabilities <- function(x, Gamma, p, delta){
  p_mat <- matrix(0, nrow=nrow(Gamma), ncol=length(x))
  for(i in 1:nrow(Gamma)){
    for(j in 1:length(x)){
      p_mat[i, j] <- p(i, x[j])
    }
  }
  return(forward_ll_cpp(log(Gamma), log(p_mat), log(delta)))
}
