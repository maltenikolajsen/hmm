#' @describeIn forward_probabilities Matrix of backward probabilities
backward_probabilities <- function(x, Gamma, p){
  p_mat <- outer(1:ncol(Gamma), x, FUN=p)
  return(backward_prob_cpp(Gamma, p_mat))
}
