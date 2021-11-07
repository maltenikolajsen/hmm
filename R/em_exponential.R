#' em.exponential
#'
#' @description A wrapper for the em function used for exponential data. Should generally not be called by the user.
#'
#' @param obs
#' @param gamma
#' @param delta
#' @param rate
#'
#' @return
#'
#' @examples
em.exponential <- function(obs, gamma, delta, rate, ...){
  m <- length(delta)
  lls <- list()
  lls_mle <- list()
  for(i in 1:m){
    lls[[i]] <- function(x, param) dexp(x, param)
    lls_mle[[i]] <- function(x, u) sum(u) / sum(x * u)
  }
  param_lls <- as.list(rate)
  out <- em(obs, gamma, delta, lls, param_lls, lls_mle, ...)
  out$parameters <- as.numeric(out$parameters)
  return(out)
}
