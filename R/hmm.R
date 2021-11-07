#' Fitting Hidden Markov Models
#'
#' @description `hmm` is used to fit (or simply define) Hidden Markov Models.
#' Most commonly, it will be used to simply fit a vector of observed emissions which come from a (known) common distribution family to a hidden Markov model, where the parameters of the underlying Markov chain as well as the emission distribution family is then fitted using the EM-algorithm.
#'
#' Alternatively, one can use custom distribution families by providing various functions (see details).
#'
#' Finally, one can also disable optimization altogether to simply obtain an 'hmm'-object which can then be assessed through various standard methods (see below).
#'
#' @param x Numeric vector of observed emissions. Can be NULL if no estimation is desired.
#' @param Gamma Initial value of transition probability matrix for underlying Markov chain.
#' @param delta Numeric vector of probabilities of the initial distribution of the Markov chain.
#' @param dist Distribution family of emissions. Can be one of 'poisson', 'normal', 'binom', 'exponential' or NULL. If NULL, user must provide own functions for densities, MLEs and random generation (see details)
#' @param ... Parameters of emission distribution as well as parameters for EM-algorithm. See details.
#' @param estimate Logical variable whether or not parameters of model (Gamma, delta and emission parameters) should be estimated by EM-algorithm. Defaults to TRUE when data is available.
#'
#' @return An object of class 'hmm'.
#' The 'hmm'-class is equipped with the following methods: `summary`, `plot`, `fitted.values`, `residuals`, `logLik` and `simulate`.
#' An object of class 'hmm' is a list containing at least the following components:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `m`          | Number of hidden states.                                          |
#' | `dist`       | Distribution family of emissions.                                 |
#' | `Gamma`      | The transition probability matrix of the underlying Markov chain. |
#' | `delta`      | The initial distribution of the underlying Markov chain.          |
#' | `parameters` | List of parameters for the emission distribution family.          |
#' | `rdists`     | List of functions used for simulating. Mostly for internal use.   |
#'
#' If `x` is not `NULL`, it will also include:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `x`          | Vector of observations.                                           |
#' | `n`          | Number of observations.                                           |
#' | `logLik`     | Log-likelihood of parameters given the observed data.             |
#' | `AIC`        | AIC of the model.                                                 |
#' | `BIC`        | BIC of the model.                                                 |
#'
#' Finally, if estimation is performed, it will also include the following:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `EMlogLik`   | Vector of log-likelihoods at each iteration in the EM-algorithm.  |
#' | `n_iter`     | Number of iterations performed in the EM-algorithm.               |
#'
#' @details
#' TODO
#'
#' @export
#'
#' @seealso TODO
#'
#' @examples TODO
hmm <- function(x, Gamma, delta, dist=NULL, ..., estimate=!is.null(x)){

  # Initialize output
  out <- list(m=length(delta),
              dist=ifelse(is.null(dist), 'Custom', dist))

  # Add observed emissions if available
  if(!is.null(x)){
    out$n <- length(x)
    out$x <- x
  }

  # First run EM-algo if desired
  if(estimate){
    fct <- paste('em', ifelse(is.null(dist), '', '.'), dist, sep='')
    em_res <- do.call(fct, list(x, Gamma, delta, ...))

    # Update according to results
    Gamma <- em_res$Gamma
    delta <- em_res$delta
    parameters <- em_res$parameters

    # Add output elements exclusive to estimation
    out$EMlogLik <- em_res$log_likelihoods
    out$n_iter <- em_res$n_iter
  }

  # If not estimated, translate options into consistent 'parameters' list
  else{
    parameters <- switch (dist,
      NULL = param_lls,
      'poisson' = lambda,
      'normal' = list(mean=mean, sd=sd),
      'binom' = list(size=size, prob=prob),
      'exponential' = rate
    )
  }

  # Calculate log-likelihood, AIC and BIC if x is not NULL
  if(!is.null(x)){
    p <- switch (dist,
     NULL = function(state, x) {lls[[state]](x, parameters[[state]])},
     'poisson' = function(state, x) {dpois(x, parameters[state])},
     'normal' = function(state, x) {dnorm(x, mean=parameters$mean[state], sd=parameters$sd[state])},
     'binom' = function(state, x) {dbinom(x, size=parameters$size, prob=parameters$prob[state])},
     'exponential' = function(state, x) {dexp(x, parameters[state])}
    )
    log_alpha <- forward_logprobabilities(x, Gamma, p, delta)
    k <- max(log_alpha[, length(x)])
    logLik <- k + log(sum(exp(log_alpha[, length(x)] - k)))

    # Calculate AIC and BIC
    n_param <- length(delta)^2-1 + length(unlist(parameters))
    AIC <- -2*logLik + 2*n_param
    BIC <- -2*logLik + n_param * log(length(x))

    # Add to output
    out$logLik <- logLik
    out$AIC <- AIC
    out$BIC <- BIC
  }

  # Make list of rdist functions
  if(!is.null(dist)){
    rdists = list()
    # NOTE: R's promise objects are *very* weird here, so some hacky fixes are required - but it works!
    rdist_state <- function(state){
      state
      switch (dist,
       'poisson' = function(n) {rpois(n, parameters[state])},
       'normal' = function(n) {rnorm(n, mean=parameters$mean[state], sd=parameters$sd[state])},
       'binom' = function(n) {rbinom(n, size=parameters$size, prob=parameters$prob[state])},
       'exponential' = function(n) {rexp(n, parameters[state])}
      )
    }

    for(j in 1:length(delta)){
      rdists[[j]] <- rdist_state(j)
    }
  }

  # Add defaults to output
  out <- c(out, list(
    Gamma=Gamma,
    delta=delta,
    parameters=parameters,
    rdists=rdists
  ))

  class(out) <- 'hmm'
  return(out)
}
