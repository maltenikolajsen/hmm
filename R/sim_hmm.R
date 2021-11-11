#' sim_hmm
#'
#' @description A function for simulating HMM's. Not expected to be called by user, use instead the 'simulate' method for hmm objects.
#'
#' @param n Integer. Length of hmm
#' @param delta m-dimensional vector. Initial distribution of Markov chain
#' @param Gamma m x m matrix. Transition probs for Markov chain
#' @param rdists List of functions f_i, where f_i(k) generates k i.i.d. copies of X|C=i
#' @param include_state Logical. Denotes whether or not the hidden states should be included. If TRUE, they are appended to the end of the HMM.
#'
#'
#' @return Vector of simulated HMM
#'
sim_hmm <- function(n, delta, Gamma, rdists, include_state=FALSE){

  # Set m for notation
  m <- ncol(Gamma)

  # First, generate C_1,...,C_n
  C <- numeric(n)
  C[1] <- sample(1:m, 1, prob=delta)

  for(i in 2:n){
    probs <- Gamma[C[i-1], ]
    C[i] <- sample(1:m, 1, prob=Gamma[C[i-1], ])
  }

  # Next, simulate X_i's
  X <- numeric(n)
  for(j in 1:m){
    X[C == j] <- rdists[[j]](sum(C == j))
  }

  if(include_state){
    return(c(X, C))
  }

  return(X)
}
