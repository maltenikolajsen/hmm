state_prediction <- function(obs,delta,gamma,lls,param_lls,h)
  {
  m <- ncol(gamma)
  n <- length(obs)

  p <- function(state, x){
    lls[[state]](x, param_lls[[state]])
  }
  p <- Vectorize(p)

  p_mat <- outer(1:m, obs, p)
  log_p_mat <- log(p_mat)

  log_delta <- log(delta)
  log_gamma <- log(gamma)
  log_alpha <- forward_logprobabilities(obs, gamma, p, delta)

  log_theta <- rep(NA, m)
  k <- max(log_alpha[, n])
  denominator <- k + log(sum(exp(log_alpha[, n] - k)))
  log_theta <- log_alpha[, n] - denominator

  log_gamma_power <- matrix(0, nrow = m, ncol = m)
  for(i in 1:h){
    log_gamma_power <- log_gamma + log_gamma_power
  }

  log_prob <- rep(NA, m)
  for(i in 1:m){
    temp <- log_gamma_power[,i] + log_theta
    k <- max(temp)
    log_prob[i] <- k + log(sum(exp(temp - k)))
  }

  exp(log_prob)
}
