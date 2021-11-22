#' Parameter estimation by EM algorithm
#'
#' @keywords internal
#'
#' @description
#' Calculates approximate MLE's of all involved parameters, i.e. initial distribution probabilities, transition probabilities and whichever parameters are required for the marginal distribution of the emissions.
#' Not expected to be called by the user directly, instead being called during the initialization of an `hmm` object.
#'
#' @param obs The observed data
#' @param gamma Initial value of transition matrix
#' @param delta Initial value of initial distribution
#' @param lls List of marginal densities
#' @param param_lls List of parameters for marginal densities
#' @param lls_mle List of functions that generate MLE's
#' @param epsilon Threshold to stop algorithm when difference in log-likelihood is less than this
#' @param max_iter Maximum number of iterations before stopping
#' @param ... Additional arguments, mainly here to avoid "unused argument"-errors.
#'
#' @return A list of log-likelihoods, parameter estimates and number of iterations.
em <- function(obs, gamma, delta, lls, param_lls, lls_mle, epsilon = 1e-5, max_iter = 1000, ...){
  ##################
  # ERROR HANDLING #
  ##################

  # Check that the number of functions matches m
  if(length(lls) != ncol(gamma) || length(lls_mle) != ncol(gamma) || length(param_lls) != ncol(gamma)){
    stop('Number of estimation functions and parameters must match to the number of states!')
  }

  ##################

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
    log_ll <- k + log(sum(exp(log_alpha[,n] - k)))

    # If log-likelihood is NaN, something has gone wrong (most likely with the provided estimates)
    # Issue warning and break loop before error occurs
    if(is.nan(log_ll)){
      warning('Log-likelihood is NaN! Aborting EM-algorithm...')
      break
    }

    # Break if change in log-likelihood is small
    if(iteration > 1 && abs(log_ll - logLs[iteration-1]) < epsilon){
      break
    }
    logLs[iteration] <- log_ll

    #E step

    # u_hat = matrix where u_hat[i, j]=P(C_j=i | X=obs)
    log_u_hat <- log_beta + log_alpha - log_ll
    u_hat <- exp(log_u_hat)

    #M step

    # Update transition probs and initial dist
    delta <- u_hat[,1]
    gamma <- exp(update_log_gamma_cpp(log_alpha, log_beta, log_gamma, log_p_mat, log_ll))

    # Update parameters for conditional densities
    for(i in 1:m){
      param_lls[[i]] <- lls_mle[[i]](obs, u_hat[i,])
    }
  }

  # Raise warning if max number of iterations reached
  if(length(logLs) == max_iter){
    warning("Maximum number of iterations reached - parameters unlikely to be MLE's")
  }

  return(list(log_likelihoods = logLs,
              n_iter = length(logLs),
              delta = delta,
              Gamma = gamma,
              parameters = param_lls))
}
