#' @title em
#'
#' @export
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


    #log_beta <- get_log_beta(log_p_mat, log_gamma)
    log_beta <- backward_ll_cpp(log_gamma, log_p_mat)
    #log_alpha <- forward_logprobabilities(p_mat, gamma, delta) // Der var fejl i disse alligevel
    log_alpha <- forward_ll_cpp(log_gamma, log_p_mat, log_delta)

    #Noget galt!
    # Ja, at der er andet end én funktion i denne R fil lol

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
