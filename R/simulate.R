#' Simulating Hidden Markov objects
#'
#' @param object Object of class `hmm`.
#' @param nsim Number of simulations.
#' @param seed Seed to be used for random generation.
#' @param include_state Logical, whether or not the hidden state should also be returned.
#' @param ... Additional arguments.
#'
#' @return If observations are available in the `hmm` object, will return a matrix with `nsim` rows, each of which being a simulation of equal length as the observed data.
#' Otherwise, will return a vector of length `nsim` of simulations from the model.
#' In either case, if `include_state` is `TRUE`, each simulation will be appended with a vector of equal length indicating the hidden state.
#'
#' @export
#'
#' @examples
#' # Continuing examples from hmm page
#' \dontshow{example(hmm)}
#' # Normal distributions
#' X.normal <- simulate(hmm.normal, nsim=100)
#' summary(X.normal)
#' plot(X.normal)
#' hist(X.normal)
#'
#' # Custom (uniform) distributions
#' Z.unif <- simulate(hmm.unif, nsim=200, include_state=TRUE)
#' X.unif <- Z.unif[1:100]
#' G.unif <- Z.unif[101:200]
#' plot(X.unif, type='h', col=G.unif, lwd=2)
#'
#' # Custom (mixed) distributions
#' Z.mixture <- simulate(hmm.mixture, nsim=200, include_state=TRUE)
#' X.mixture <- Z.mixture[1:200]
#' G.mixture <- Z.mixture[201:400]
#' plot(X.mixture, type='h', col=G.mixture, lwd=2)
simulate.hmm <- function(object, nsim = 1, seed = NULL, include_state = FALSE, ...){
  if(!is.null(seed)){set.seed(seed)}

  # Evil hack because '...' is weird
  dots <- list(...)
  include_state <- dots$include_state
  if(is.null(include_state)){include_state <- FALSE}

  # Seperate into cases where observations are available and not
  if(is.null(object$x)){
    return(sim_hmm(nsim, object$delta, object$Gamma, object$rdists, include_state))
  }

  out <- matrix(NA, nrow = nsim, ncol = ifelse(include_state, 2*object$n, object$n))
  for(j in 1:nsim){
    out[j,] <- sim_hmm(object$n, object$delta, object$Gamma, object$rdists, include_state)
  }
  return(out)
}
