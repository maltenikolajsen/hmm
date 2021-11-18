#' Global state prediction using Viterbi's algorithm
#'
#' @keywords internal
#'
#' @description A function for global decoding of HMM's. Not expected to be called by user. Called by the general HMM function.
#'
#' @param obs Observations.
#' @param delta m-dimensional vector. Initial distribution of Markov chain
#' @param Gamma m x m matrix. Transition probs for Markov chain
#' @param p Function returning the evaluated density given the state.
#'
#'
#' @return List of most likely states and decoding probabilities.
#'
viterbi <- function(obs, delta, Gamma, p){

  n <- length(obs); m <- ncol(Gamma)
  xi <- matrix(NA, ncol = n, nrow = m)
  foo <- delta * sapply(1:m, function(k) p(k, obs[1]))

  xi[,1] <- foo / sum(foo)
  for(i in 2:n){
    foo <- apply(xi[,i-1]*Gamma,2,max) * sapply(1:m, function(k) p(k, obs[i]))
    xi[,i] <- foo/sum(foo)
  }


  iv <- numeric(n)
  probs <- matrix(NA, nrow = m, ncol = n)
  probs[,n] <- xi[,n]
  iv[n] <- which.max(probs[,n])
  for(i in (n-1):1){
    probs[,i] <- Gamma[,iv[i+1]]*xi[,i]
    iv[i] <- which.max(probs[,i])
  }
  return(list(states = iv, probs = probs))
}

