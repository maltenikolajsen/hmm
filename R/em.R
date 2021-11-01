get_log_beta <- function(
  log_p_mat,
  log_gamma
){
  m <- ncol(log_gamma); n <- length(log_p_mat) / m

  log_beta <- matrix(NA, nrow = m, ncol = n)

  log_beta[,n] <- rep(0, m)

  for(i in (n-1):1){
    tmp_logbeta <- rep(NA, m)
    for(j in 1:m){
      #tmp_zeta <- log(transition[j,]) + log(helper_function(observation[i])) + logbeta[[i+1]]

      tmp_zeta <- log_gamma[j,] + log_p_mat[,i] + log_beta[,i+1]
      logxi_star <- max(tmp_zeta)
      tmp_logbeta[j] <- logxi_star + log(sum(exp(tmp_zeta - logxi_star)))
    }
    log_beta[,i] <- tmp_logbeta
  }

  log_beta
}

get_log_alpha <- function(
  log_p_mat,
  log_gamma,
  log_delta
){
  m <- ncol(log_gamma); n <- length(log_p_mat) / m

  log_alpha <- matrix(NA, nrow = m, ncol = n)

  log_alpha[,1] <-  log_p_mat[,1] + log_delta

  for(i in 2:n){
    for(j in 1:m){
      tmp_zeta <- as.vector(log_gamma[j,]) + log_p_mat[,i] + log_alpha[,i-1]
      logxi_star <- max(tmp_zeta)
      log_alpha[j,i] <- logxi_star + log(sum(exp(tmp_zeta - logxi_star)))
    }
  }

  log_alpha
}


em <- function(obs, gamma, delta, lls, param_lls, lls_mle, max_iter = 1){

  n <- length(obs)
  m <- ncol(gamma)

  p <- function(state, x){
    lls[[state]](x, param_lls[[state]])
  }
  p <- Vectorize(p)

  for(iteration in 1:max_iter){

    p_mat <- outer(1:m, obs, p)
    log_p_mat <- log(p_mat)
    log_gamma <- log(gamma)
    log_delta <- log(delta)


    log_beta <- get_log_beta(log_p_mat, log_gamma)
      #backward_logprobabilities(p_mat, gamma)
    log_alpha <- forward_logprobabilities(p_mat, gamma, delta)
      #forward_logprobabilities(p_mat, gamma, delta)

    #Noget galt!

    c <- max(log_alpha[,n])
    log_ll <- c + log(sum(exp(log_alpha[,n] - c)))

    #E step
    log_u_hat <- log_beta + log_alpha - log_ll
    u_hat <- exp(log_u_hat)

    log_f <- matrix(NA, nrow = m, ncol = m)
    f <- matrix(NA, nrow = m, ncol = m)
    for(j in 1:m){
      for(k in 1:m){
        log_v_hat_jk <- log_alpha[j,1:(n-1)] + log_gamma[j,k] + log_p_mat[k,2:n] + log_beta[k,2:n] - log_ll
        v_hat_jk <- exp(log_v_hat_jk)
        #print("NEXT")
        #print(log_v_hat_jk)
        #print(log_v_hat_jk)
        #print("NEXT")
        f[j,k] <- sum(v_hat_jk)
        c <- max(log_v_hat_jk)
        log_f[j,k] <- c + log(sum(exp(log_v_hat_jk - c)))
      }
    }


    #M step
    delta <- u_hat[,1] / sum(u_hat[,1])
    gamma <- f / rowSums(f)

    for(i in 1:m){
      param_lls[[i]] <- lls_mle[[i]](obs, u_hat[i,])
    }

    #print(param_lls)
    #print(-log_ll)
  }
  list(delta = delta, gamma = gamma, param_lls = param_lls, log_ll = log_ll)
  #list(log_alpha, p_mat)
  #log_alpha
}


obs <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
delta <- c(.3, .7)
gamma <- matrix(c(.25, .75, .75, .25), nrow = 2, byrow = T)
lls <- list(function(x, param) dpois(x, param), function(x, param) dpois(x, param))
param_lls <- list(10, 30)
lls_mle <- list(function(x, u) sum(x * u) / sum(u), function(x, u) sum(x * u) / sum(u))
a <- em(obs, gamma, delta, lls, param_lls, lls_mle, max_iter = 20)
print(a)
