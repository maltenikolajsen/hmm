#' @title em
#'
#' @export
em <- function(obs, gamma, delta, lls, param_lls, lls_mle, max_iter = 1){

  # Set n and m for notation
  n <- length(obs)
  m <- ncol(gamma)

  # Define p-function in proper format for backward- and forward algos
  p <- function(state, x){
    lls[[state]](x, param_lls[[state]])
  }
  p <- Vectorize(p)

  # Iterate E- and M-step a total of max_iter times
  # TODO: Add convergence break
  for(iteration in 1:max_iter){

    # Get forward and backward log-probabilities
    p_mat <- outer(1:m, obs, p)
    log_p_mat <- log(p_mat)
    log_gamma <- log(gamma)
    log_delta <- log(delta)

    log_beta <- backward_ll_cpp(log_gamma, log_p_mat)
    log_alpha <- forward_ll_cpp(log_gamma, log_p_mat, log_delta)

    # Get log-likelihood of entire data-set
    k <- max(log_alpha[,n])
    log_ll <- k + log(sum(exp(log_alpha[,n] - k))) # Please don't use 'c' as variable

    #E step

    # u_hat = matrix where u_hat[i, j]=P(C_j=i | X=obs)
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
    delta <- u_hat[,1]
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
