#' em.normal
#'
#' @description A wrapper for the em function used for normal data. Should generally not be called by the user.
#'
#' @param obs
#' @param gamma
#' @param delta
#' @param mu
#' @param sigma
#'
#' @return
#'
#' @examples
em.normal <- function(obs, gamma, delta, mu, sigma, ...){
  m <- length(delta)
  lls <- list()
  lls_mle <- list()
  param_lls <- list()
  for(i in 1:m){
    lls[[i]] <- function(x, param) dnorm(x, mean=param[1], sd=param[2])
    lls_mle[[i]] <- function(x, u) {
      mu_hat <- sum(u*x) / sum(u)
      sigma_hat <- sqrt(sum(u*(x-mu_hat)^2) / sum(u))
      return(c(mu_hat, sigma_hat))
    }
    param_lls[[i]] <- c(mu[i], sigma[i])
  }
  return(em(obs, gamma, delta, lls, param_lls, lls_mle, ...))
}
