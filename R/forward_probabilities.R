#' forward_probabilities
#'
#' @description Calculate matrix of forward probabilities for a given data set and parameters in a HMM. Essentially just a wrapper for the true function written in C++.
#'
#' @param x Vector of observed emissions
#' @param Gamma Matrix (m x m) of transition probabilities of the underlying Markov chain
#' @param p A function where p(i, x)=P(X_t=x | Y_t=i)
#' @param delta Vector of starting probabilities, i.e. the distribution of Y_1
#'
#' @details Let (X, Y)_i, i=1,...,n constitute a HMM, that is, Y_1,...,Y_n is a Markov chain with states 1,...,m and the X_i's are dependent only through the Y_i's with common conditional distributions X_i | Y_i=j = P_j for all i=1,...,n. This function then generates a matrix where the i,j'th entry is P(X_1=x_1, ..., X_j=x_j, Y_j=i), where x_1,...,x_n are the observed values (emissions) of X.
#'
#' @return A matrix (m x n) of forward probabilities.
#' @export
#'
#' @examples TODO
forward_probabilities <- function(x, Gamma, p, delta){
  p_mat <- outer(1:ncol(Gamma), x, FUN=p)
  return(forward_prob_cpp(Gamma, p_mat, delta))
}
