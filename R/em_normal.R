#' em.normal
#'
#' @description A wrapper for the em function used for normal data. Should generally not be called by the user.
#'
#' @param obs
#' @param gamma
#' @param delta
#' @param mean
#' @param sd
#'
#' @return
#'
#' @examples
em.normal <- function(obs, gamma, delta, mean, sd, ...){
  m <- length(delta)
  lls <- list()
  lls_mle <- list()
  param_lls <- list()
  for(i in 1:m){
    lls[[i]] <- function(x, param) dnorm(x, mean=param[1], sd=param[2])
    lls_mle[[i]] <- function(x, u) {
      mean_hat <- sum(u*x) / sum(u)
      sd_hat <- sqrt(sum(u*(x-mean_hat)^2) / sum(u))
      return(c(mean_hat, sd_hat))
    }
    param_lls[[i]] <- c(mean[i], sd[i])
  }
  out <- em(obs, gamma, delta, lls, param_lls, lls_mle, ...)
  out$parameters <- list(mean=unlist(out$parameters)[c(TRUE, FALSE)], # Weird hack, but it works
                         sd=unlist(out$parameters)[c(FALSE, TRUE)])
  return(out)
}
