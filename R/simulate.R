#' Simulate method for class 'hmm'.
#' @export
simulate.hmm <- function(model, nsim = 1, seed = NULL, ...){
  set.seed(seed)

  # Evil hack because '...' is weird
  dots <- list(...)
  include_state <- dots$include_state
  if(is.null(include_state)){include_state <- FALSE}

  # Seperate into cases where observations are available and not
  if(is.null(model$x)){
    return(sim_hmm(nsim, model$delta, model$Gamma, model$rdists, include_state))
  }

  out <- matrix(NA, nrow = nsim, ncol = ifelse(include_state, 2*model$n, model$n))
  for(j in 1:nsim){
    out[j,] <- sim_hmm(model$n, model$delta, model$Gamma, model$rdists, include_state)
  }
  return(out)
}
