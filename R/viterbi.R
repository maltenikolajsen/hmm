#' viterbi
#'
#' @description A function for global decoding of HMM's. Not expected to be called by user. Called by the general HMM function.
#'
#' @param obs Observations.
#' @param delta m-dimensional vector. Initial distribution of Markov chain
#' @param Gamma m x m matrix. Transition probs for Markov chain
#' @param p Function returning the evaluated density given the state.
#'
#'
#' @return Vector of hidden states most likely observed using global decoding.
#'
viterbi <- function(obs, delta, Gamma, p){

  n <- length(obs); m <- ncol(Gamma)
  xi <- matrix(NA, ncol = n, nrow = m)
  foo <- delta * p(1:m, obs[1])

  xi[,1] <- foo / sum(foo)
  for(i in 2:n){
    foo <- apply(xi[,i-1]*Gamma,2,max) * p(1:m, obs[i])
    xi[,i] <- foo/sum(foo)
  }

  iv <- numeric(n)
  iv[n] <- which.max(xi[,n])
  for(i in (n-1):1){
    iv[i] <- which.max(Gamma[,iv[i+1]]*xi[,i])
  }
  iv
}

