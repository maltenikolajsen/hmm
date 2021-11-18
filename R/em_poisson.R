#' @describeIn em EM-algorithm for Poisson emissions
em.poisson <- function(obs, gamma, delta, lambda, ...){
  m <- length(delta)
  lls <- list()
  lls_mle <- list()
  for(i in 1:m){
    lls[[i]] <- function(x, param) dpois(x, param)
    lls_mle[[i]] <- function(x, u) sum(x * u) / sum(u)
  }
  param_lls <- as.list(lambda)
  out <- em(obs, gamma, delta, lls, param_lls, lls_mle, ...)
  out$parameters <- as.numeric(out$parameters)
  return(out)
}
