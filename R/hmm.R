#' Fitting Hidden Markov Models
#'
#' @description `hmm` is used to fit (or simply define) Hidden Markov Models.
#' Most commonly, it will be used to simply fit a vector of observed emissions which come from a (known) common distribution family to a hidden Markov model, where the parameters of the underlying Markov chain as well as the emission distribution family is then fitted using the EM-algorithm.
#'
#' Alternatively, one can use custom distribution families by providing various functions (see details).
#'
#' Finally, one can also disable optimization altogether to simply obtain an 'hmm'-object which can then be assessed through various standard methods (see below).
#'
#' @param x Numeric vector of observed emissions. Can be `NULL` if no estimation is desired.
#' @param Gamma Initial value of transition probability matrix for underlying Markov chain.
#' @param delta Numeric vector of probabilities of the initial distribution of the Markov chain.
#' @param dist Distribution family of emissions. Can be one of 'poisson', 'normal', 'binom', 'exponential' or 'custom'. If 'custom', user must provide own functions for densities, MLEs and random generation (see details)
#' @param ... Parameters of emission distribution as well as parameters for EM-algorithm. See details.
#' @param estimate Logical variable whether or not parameters of model (Gamma, delta and emission parameters) should be estimated by EM-algorithm. Defaults to `TRUE` when data is available.
#'
#' @return An object of class 'hmm'.
#' The 'hmm'-class is equipped with a variety of default methods, see 'See also' section for details.
#' An object of class 'hmm' is a list containing at least the following components:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `m`          | Number of hidden states.                                          |
#' | `dist`       | Distribution family of emissions.                                 |
#' | `Gamma`      | The transition probability matrix of the underlying Markov chain. |
#' | `delta`      | The initial distribution of the underlying Markov chain.          |
#' | `parameters` | List of parameters for the emission distribution family.          |
#' | `rdists`     | List of functions used for simulating. Only for internal use.     |
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
#' | `viterbi_s`  | States of global decoding of the model.                           |
#' | `posterior_s`| States of local decoding of the model.                            |
#' | `viterbi_s`  | Probs. of global decoding of the model.                           |
#' | `posterior_s`| Probs. of local decoding of the model.                            |
#'
#' Finally, if estimation is performed, it will also include the following:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `EMlogLik`   | Vector of log-likelihoods at each iteration in the EM-algorithm.  |
#' | `n_iter`     | Number of iterations performed in the EM-algorithm.               |
#'
#' @details
#' ## Theory
#'
#' This function is used to fit or define a hidden Markov model, i.e. a model where the distributions of X_1, ..., X_n (the
#' emissions) depend on a hidden sequence Y_1, ..., Y_n (the hidden states).
#' In particular, in this model, Y_1, ..., Y_n constitutes a Markov chain on the finite state space {1, ..., m}, and the
#' distribution of X_i depends only on the value of Y_i, i.e. X_i | Y_i=k ~ X_j | Y_j=k for all i, j.
#'
#' As such, to estimate parameters in this model, one must estimate not only the parameters of the m emission distributions,
#' but also the parameters of the Markov chain, i.e. the initial distribution delta (that is, Y_1 ~ delta) and the transition
#' matrix Gamma (that is, P(Y_k=j | Y_(k-1)=i)=Gamma_i,j)
#'
#' ## Argument details
#'
#' The `dist` and `...` arguments determine the distribution of the emissions in the different states.
#' If `dist` is one of 'poisson', 'normal', 'binom' or 'exponential', distribution parameters must be passed to the `...` argument.
#' In particular they must have the following format:
#'
#' |               |                |                                                                                                                                                      |
#' |---------------|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
#' | 'poisson'     | `lambda`       | A vector of Poisson parameters, one per hidden state.                                                                                                |
#' | 'normal'      | `mean`, `sd`   | Vectors of means and standard deviations, one per hidden state.                                                                                      |
#' | 'binom'       | `size`, `prob` | `size` is a single integer, denoting the common size for all hidden states, while `prob` is a vector of success probabilities, one per hidden state. |
#' | 'exponential' | `rate`         | Vector of rates for the exponential distributions, one per hidden state.                                                                             |
#'
#' If `dist` is 'custom', user must provide the following:
#'
#' |              |                                                                                                                                                                                                                                                                                                                                  |
#' |--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
#' | `lls`        | Function or list of functions, where each function is on the form `f(x, param)`, and returns the conditional density of X_i given the parameter `param`.                                                                                                                                                                         |
#' | `params_lls` | List of parameters, one for each state. Must fit with `lls` in the sense that `lls(x, param[[i]])` (or `lls[[i]](x, param[[i]])` if `lls` is a list) denotes the density in point x of X_j given Y_j=i                                                                                                                           |
#' | `lls_mle`    | Function or list of functions, which return the maximum likelihood estimates given the provided data `x` and a vector (of equal length) of scalars `u` where each scalar is between 0 and 1. That is, the function(s) must be on the form `h(x, u)` and must return the value of `param` maximizing `sum(u * log(f(x, param)))`. |
#' | `rdist`      | Function or list of functions that generate random emissions from the different hidden states. The function(s) must be on the form `r(n, param)` and must return a vector of `n` random realizations given the parameter `param`.                                                                                                |
#'
#' For each of `lls`, `lls_mle` and `rdist`, if the distribution family itself (and not just the parameters) depends on the hidden state, a list of functions must be provided, one for each state.
#'
#' Finally, `...` also takes arguments passed to the EM-algorithm, namely: `epsilon`, .
#'
#' |            |                                                                                                         |
#' |------------|---------------------------------------------------------------------------------------------------------|
#' | `epsilon`  | The value at which estimation halts, if the difference in log-likelihood between iterations falls below |
#' | `max_iter` | Maximum number of iterations before halting.                                                            |
#'
#'
#' @export
#'
#' @seealso The `hmm` object has methods for the following generic functions:
#' [AIC][AIC.hmm], [BIC][BIC.hmm], [fitted.values][fitted.hmm], [logLik][logLik.hmm], [plot][plot.hmm], [print][print.hmm], [residuals][residuals.hmm], [simulate][simulate.hmm] and [summary][summary.hmm].
#' Some (most) of these are only available, when data is provided, i.e. when `x` is not `NULL`.
#'
#' @examples
#' # Annual counts of earthquakes magnitude 7 or greater, 1900-2006.
#' # Source:
#' # Earthquake Data Base System of the U.S. Geological Survey, National
#' # Earthquake Information Center, Golden CO
#'
#' quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
#' Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
#' delta <- c(1, 1)/2
#' lambda <- c(10, 30)
#' hmm.EQ <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda)
#' hmm.EQ
#'
#' # If one does not want estimation by EM algorithm (e.g. for comparison of summary statistics), it can be disabled
#' hmm.EQ_no_opt <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda, estimate=FALSE)
#' hmm.EQ_no_opt
#'
#' # Creating 'empty' hmm object for sake of simulation (see simulate for further details)
#' # Here where all emission distributions are normal
#' Gamma <- rbind(c(0.5, 0.25, 0.25),
#'                c(0.1, 0.8 , 0.1),
#'                c(  0, 0.2 , 0.8))
#' delta <- c(1, 0, 0)
#' mean <- c(0, 5, 10)
#' sd <- rep(1, 3)
#'
#' hmm.normal <- hmm(NULL, Gamma=Gamma, delta=delta, dist='normal', mean=mean, sd=sd)
#' hmm.normal
#'
#' # Here, the emission distributions are custom (Uniform[0, theta])
#' Gamma <- rbind(c(0.5, 0.25, 0.25),
#'                c(0.1, 0.8 , 0.1),
#'                c(  0, 0.2 , 0.8))
#' delta <- c(1, 0, 0)
#' theta <- list(1, 5, 10)
#' lls <- function(x, param){dunif(x, 0, param)}
#' lls_mle <- function(x, u){max(x)}
#' rdist <- function(n, param){runif(n, 0, param)}
#' hmm.unif <- hmm(NULL, Gamma=Gamma, delta=delta, lls=lls, param_lls=theta, lls_mle=lls_mle, rdist=rdist)
#' hmm.unif
#'
#' # Here, the emission distributions is either normal(0, 1) or exponential(1)
#' Gamma <- rbind(c(0.2, 0.8),
#'                c(0.8, 0.2))
#' delta <- c(1, 1)/2
#' param <- list(c(0, 1), 1)
#'
#' lls <- list(function(x, param){dnorm(x, param[1], param[2])},
#'             function(x, param){dexp(x, param)})
#'
#' lls_mle <- list(function(x, u){mean_hat <- sum(u*x) / sum(u); c(mean_hat, sqrt(sum(u*(x-mean_hat)^2) / sum(u)))},
#'                 function(x, u){sum(u)/sum(u*x)})
#'
#' rdist <- list(function(n, param){do.call(rnorm, args=as.list(c(n, param)))},
#'               function(n, param){rexp(n, param)})
#'
#' hmm.mixture <- hmm(NULL, Gamma=Gamma, delta=delta, lls=lls, param_lls=param, lls_mle=lls_mle, rdist=rdist)
#' hmm.mixture
hmm <- function(x, Gamma, delta, dist='custom', ..., estimate=!is.null(x)){
  ##################
  # ERROR HANDLING #
  ##################

  # Check that Gamma and delta match in length/dimension and that they are valid
  if(ncol(Gamma) != nrow(Gamma) || !all((rowSums(Gamma) - 1) < 1e-5)){
    stop('Gamma is not a valid transition matrix! (It must be square and have row sums = 1)')
  }

  if(length(delta) != ncol(Gamma) || sum(delta) != 1){
    stop('delta is not a valid probability vector! (It must be of length m and sum to 1)')
  }

  # Check that we don't estimate if x=NULL or length(x)<2
  if((is.null(x) || length(x) < 2) && estimate){
    stop('Need at least 2 data points to estimate!')
  }

  # Check that dist is valid
  if(!dist %in% c('custom', 'poisson', 'normal', 'binom', 'exponential')){
    stop('Invalid distribution! Must be one of "custom", "poisson", "normal", "binom" or "exponential"')
  }

  ##################

  # Initialize output
  out <- list(m=length(delta),
              dist=dist)

  # If custom dist is provided, but only one family, make lists of appropriate parameters
  if(dist == 'custom'){
    if(typeof(lls) == 'closure'){lls <- rep(list(lls), length(delta))}
    if(typeof(lls_mle) == 'closure'){lls_mle <- rep(list(lls_mle), length(delta))}
    if(typeof(rdist) == 'closure'){rdist <- rep(list(rdist), length(delta))}
  }

  # Add observed emissions if available
  if(!is.null(x)){
    out$n <- length(x)
    out$x <- x
  }

  # First run EM-algo if desired
  if(estimate){
    fct <- ifelse(dist=='custom', 'em', paste('em.', dist, sep=''))
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
    # There is a *very* strange behaviour in R where only the last argument in
    # '...' is available in the function's environment (even though hasArg
    # returns TRUE????), so it must be unpacked in a list.
    params.list <- list(...)
    parameters <- switch (dist,
      'custom' = params.list$param_lls,
      'poisson' = params.list$lambda,
      'normal' = list(mean=params.list$mean, sd=params.list$sd),
      'binom' = list(size=params.list$size, prob=params.list$prob),
      'exponential' = params.list$rate
    )
  }

  # Calculate log-likelihood, AIC and BIC if x is not NULL
  if(!is.null(x)){
    p <- switch (dist,
     'custom' = function(state, x) {lls[[state]](x, parameters[[state]])},
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

    # Add global decoding if possible.
    if(!is.nan(logLik)){
      global_decoding <- viterbi(obs = x, delta = delta, Gamma = Gamma, p = p)
      local_decoding <- l_decoding(obs =  x, delta = delta, Gamma = Gamma, p = p)

      out$viterbi_s <- global_decoding$states
      out$posterior_s <- local_decoding$states

      out$viterbi_p <- global_decoding$probs
      out$posterior_p <- local_decoding$probs
    }
  }

  # Make list of rdist functions
  rdists = list()
  # NOTE: R's promise objects are *very* weird here, so some hacky fixes are required - but it works!
  rdist_state <- function(state){
    state
    switch (dist,
     'custom' = function(n) {rdist[[state]](n, parameters[[state]])},
     'poisson' = function(n) {rpois(n, parameters[state])},
     'normal' = function(n) {rnorm(n, mean=parameters$mean[state], sd=parameters$sd[state])},
     'binom' = function(n) {rbinom(n, size=parameters$size, prob=parameters$prob[state])},
     'exponential' = function(n) {rexp(n, parameters[state])}
    )
  }

  for(j in 1:length(delta)){
    rdists[[j]] <- rdist_state(j)
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
