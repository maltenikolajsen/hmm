#' em
#'
#' @description
#' Calculates approximate MLE's of all involved parameters, i.e. initial distribution probabilities, transition probabilities and whichever parameters are required for the marginal distribution of the emissions.
#'
#' @param obs The observed data
#' @param gamma Initial value of transition matrix
#' @param delta Initial value of initial distribution
#' @param lls List of marginal densities
#' @param param_lls List of parameters for marginal densities
#' @param lls_mle List of functions that generate MLE's
#' @param epsilon Threshold to stop algorithm when difference in log-likelihood is less than this
#' @param max_iter Maximum number of iterations before stopping
#'
#' @details
#' Add some details.
#'
#' @return A list of log-likelihoods, parameter estimates and number of iterations.
#' @export
#'
#' @examples
#' #TODO
#'
em <- function(obs, gamma, delta, lls, param_lls, lls_mle, epsilon = 10^(-4), max_iter = 1000){

  # Create vector of log-likelihoods - initialize with Inf for comparison purposes
  logLs <- c()

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

    # Break if change in log-likelihood is small
    if(iteration > 1 && abs(log_ll - logLs[iteration-1]) < epsilon){
      break
    }
    logLs[iteration] <- log_ll

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
        k <- max(log_v_hat_jk)
        log_f[j,k] <- k + log(sum(exp(log_v_hat_jk - k)))
      }
    }


    #M step

    # Update transition probs and initial dist
    delta <- u_hat[,1]
    gamma <- f / rowSums(f)

    # Update parameters for conditional densities
    for(i in 1:m){
      param_lls[[i]] <- lls_mle[[i]](obs, u_hat[i,])
    }
  }
  return(list(log_likelihoods = logLs,
              n_iter = length(logLs),
              delta = delta,
              gamma = gamma,
              param_lls = param_lls,
              log_likelihood = log_ll))
}
