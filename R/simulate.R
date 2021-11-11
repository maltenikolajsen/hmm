#' Simulate method for class 'hmm'.
#' @export
simulate.hmm <- function(model, nsim = 1, seed = NULL, ...){
  set.seed(NULL)
  hasArg(include_state) || (include_state = F)

  out <- matrix(NA, nrow = nsim, ncol = ifelse(include_state, model$n, 2*model$n))
  for(j in 1:nsim){
    out[j,] <- sim_hmm(model$n, model$delta, model$Gamma, model$rdists, include_state = F)
  }
  return(out)
}
