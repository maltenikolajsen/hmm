#' em.poisson
#'
#' @description A wrapper for the em function used for Poisson data. Should generally not be called by the user.
#'
#' @param obs
#' @param gamma
#' @param delta
#' @param lambda
#'
#' @return
#'
#' @examples
em.poisson <- function(obs, gamma, delta, lambda, ...){
  m <- length(delta)
  lls <- list()
  lls_mle <- list()
  for(i in 1:m){
    lls[[i]] <- function(x, param) dpois(x, param)
    lls_mle[[i]] <- function(x, u) sum(x * u) / sum(u)
  }
  param_lls <- as.list(lambda)
  return(em(obs, gamma, delta, lls, param_lls, lls_mle, ...))
}
