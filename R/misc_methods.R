#' Parameter table
#'
#' @description Helper function for other methods that prints parameters of `hmm` class in a human-readable format.
#'
#' @param x Object of type `hmm`.
#'
#' @return NULL
param_table <- function(x){
  switch (x$dist,
          'custom' = {
            nparm <- max(unlist(lapply(x$parameters, length)))
            P <- c()
            for(p in x$parameters){
              P <- rbind(P, c(p, rep(NA, nparm-length(p))))
            }
            colnames(P) <- sapply(1:nparm, function(i) {paste('Param.', i)})
          },
          'poisson' = {
            P <- as.matrix(x$parameters)
            colnames(P) <- 'lambda'
          },
          'normal' = {
            P <- cbind(x$parameters$mean, x$parameters$sd)
            colnames(P) <- c('mean', 'sd')
          },
          'binom' = {
            P <- as.matrix(x$parameters$prob)
            colnames(P) <- 'prob'
          },
          'exponential' = {
            P <- as.matrix(x$parameters)
            colnames(P) <- 'rate'
          }
  )
  rownames(P) <- 1:x$m
  print(P, na.print='')
}

#' Print `hmm`
#'
#' @description Print method for class `hmm`.
#' Prints a short summary of the involved parameters, namely the initial distribution, transition matrix and emission parameters.
#'
#' @param x Object of class `hmm`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
print.hmm <- function(x, ...){
  cat('Hidden Markov Model with', x$m, 'states.\n')
  cat('Initial distribution of Markov chain is:\n')
  print(zapsmall(x$delta))
  cat('\nTransition probability matrix of Markov chain is:\n')
  colnames(x$Gamma) <- 1:x$m
  rownames(x$Gamma) <- 1:x$m
  print(x$Gamma)
  cat('\nEmission distribution family is', ifelse(x$dist == 'binom',
                                                  paste('binomial with common size', x$parameters$size, 'and probabilities:\n'),
                                                  paste(x$dist, 'with parameters:\n')))
  param_table(x)
  invisible(x)
}

#' Evaulating Hidden Markov Models
#'
#' Functions for evaluating model fitness of a hidden Markov model using AIC, BIC and log-likelihood.
#'
#' The log-likelihood is calculated using a combination of the forward algorithm for finding the (regular) likelihood and the log-sum-exp-trick to properly convert this to a log-likelihood in order to avoid underflow.
#' The AIC and BIC are calculated as usual, where the number of parameters is (m^2-1) + k where m is the number of states and k is the number of parameters in the emission distributions.
#' This is because we estimate m*(m-1) probabilities in the transition matrix and m-1 probabilities in the initial distribution vector, i.e. (m+1)(m-1) = m^2-1 in all.
#' All of this, however, is done during the creation of the `hmm` instance, so this function simply returns the values stored in the object.
#'
#' @export
#'
#' @param object Object of class `hmm`.
#' @param ... Optinally more fitted model objects.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC
#'
#' @describeIn AIC.hmm Returns the AIC of the HMM.
#'
#' @examples
#' # Continuation of Earthquake data example
#' quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
#' Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
#' delta <- c(1, 1)/2
#' lambda <- c(10, 30)
#' M <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda)
#'
#' AIC(M)
#' BIC(M)
#' logLik(M)
#'
AIC.hmm <- function(object, ..., k = 2){
  if(is.null(object$x)){
    stop('Must have at least one observation to produce AIC!')
  }
  n_param <- object$m^2-1 + length(unlist(object$parameters))
  return(-2*object$logLik + k*n_param)
}

#' @describeIn AIC.hmm Returns the BIC of the HMM.
#' @export
BIC.hmm <- function(object, ...){
  if(is.null(object$x)){
    stop('Must have at least one observation to produce BIC!')
  }
  return(object$BIC)
}

#' @describeIn AIC.hmm Returns the log-likelihood of the HMM
#' @export
logLik.hmm <- function(object, ...){
  if(is.null(object$x)){
    stop('Must have at least one observation to produce log-likelihood!')
  }
  return(object$logLik)
}

#' Fitted values for Hidden Markov Models
#'
#' @description Outputs predicted values of Hidden Markov Model in the form of a list where the first element is a vector of state predictions and second element is a vector of emission predictions (expected value corresponding predicted state).
#'
#' @param x Object of type `hmm`.
#' @param method Method for state prediction. One of either 'local' (default, resulting in local decoding), 'viterbi' or 'global' (both resulting in global decoding using Viterbi's algorithm).
#' @param ... Other arguments.
#'
#' @return List of predicted states and emission values.
#'
#' @details If the distribution family is one of the included standard families (Poisson, normal, binomial or exponential), the expected value will be the plug-in estimate using the fitted parameters of each state.
#' Otherwise, the empirical mean will be used as the expected value.
#'
#' @export
#'
#' @examples
#' # Continuation of Earthquake data example
#' quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
#' Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
#' delta <- c(1, 1)/2
#' lambda <- c(10, 30)
#' M <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda)
#'
#' fitted.values(M)
fitted.hmm <- function(x, method='local', ...){
  # Check that data is available
  if(is.null(x$x)){
    stop('Cannot estimate when no data is available!')
  }

  # Find predicted states according to method chosen
  if(method=='local'){states <- x$posterior_s}
  else{states <- x$viterbi_s}

  # Get expected values for each state
  expected <- c()
  for(i in 1:x$m){
    expected[i] <- switch (x$dist,
      'custom' = mean(x$x[states == i]),
      'poisson' = x$parameters[i],
      'normal' = x$parameters$mean[i],
      'binom' = x$parameters$size * x$parameters$prob[i],
      'exponential' = 1/x$parameters[i]
    )
  }

  # Set expected emissions in a vector
  expected_emission <- numeric(x$n)
  for(i in 1:x$m){
    expected_emission[states == i] <- expected[i]
  }

  return(list(state=states,
              emission=expected_emission))
}

#' Residuals of Hidden Markov Models
#'
#' @description Returns residuals of the observed emissions.
#'
#' @param object Object of type `hmm`.
#' @param ... Other arguments passed to \link[hmm]{fitted.hmm}.
#'
#' @return Vector of residuals of the observed emissions.
#'
#' @details If the distribution family is one of the included standard families (Poisson, normal, binomial or exponential), the expected value will be the plug-in estimate using the fitted parameters of each state.
#' Otherwise, the empirical mean will be used as the expected value.
#'
#' @export
#'
#' @examples
residuals.hmm <- function(object, ...){
  return(object$x - fitted(object)$emission)
}

#' Summarizing Hidden Markov Models
#'
#' @description Summary method for class "`hmm`"
#'
#' @param object An object of class "`hmm`", usually a result of call to \link[hmm]{hmm}.
#' @param ... Additional arguments passed to \link[hmm]{fitted.hmm}
#'
#' @details
#' The function prints a summary of the parameters of the given model in a more human-readable format.
#' If fitted, it also provides model fitness statistics in the form of log-likelihood, AIC and BIC.
#'
#' @return If no data is available, essentially returns just the object ordered in a more reasonable manner and without internal functions.
#' Otherwise, also returns fitted values and residuals.
#'
#' @export
#'
#' @examples TBD
summary.hmm <- function(object, ...){
  if(is.null(object$x)){
    out <- list(dist=object$dist,
                m=object$m,
                delta=object$delta,
                Gamma=object$Gamma,
                parameters=object$parameters)
  }
  else{
    out <- list(n=object$n,
                x=object$x,
                dist=object$dist,
                m=object$m,
                delta=object$delta,
                Gamma=object$Gamma,
                parameters=object$parameters,
                fitted=fitted(object, ...),
                resid=residuals(object),
                logLik=object$logLik,
                AIC=object$AIC,
                BIC=object$BIC)
    if(!is.null(object$n_iter)){
      out <- c(out, list(n_iter=object$n_iter))
    }
  }
  class(out) <- 'summary.hmm'
  return(out)
}

#' @inherit summary.hmm
#' @export
#'
print.summary.hmm <- function(x, ...){
  # If no data is available, just print the same as normal
  if(is.null(x$x)){
    print.hmm(x, ...)
  }

  # Otherwise, include some other shit
  cat('Hidden Markov Model with', x$m, 'states.\n')
  cat('Residuals:\n')
  if(length(x$resid) < 5){
    print(x$resid)
  }
  else{
    print(quantile(x$resid))
  }
  cat('\nInitial distribution of Markov chain is:\n')
  print(zapsmall(x$delta))
  cat('\nTransition probability matrix of Markov chain is:\n')
  colnames(x$Gamma) <- 1:x$m
  rownames(x$Gamma) <- 1:x$m
  print(x$Gamma)
  cat('\nEmission distribution family is', ifelse(x$dist == 'binom',
                                                  paste('binomial with common size', x$parameters$size, 'and probabilities:\n'),
                                                  paste(x$dist, 'with parameters:\n')))
  param_table(x)
  if(!is.null(x$n_iter)){
    cat('\nModel was fitted using the EM-algorithm in', x$n_iter, 'steps.')
  }
  cat('\nThe log-likelihood of the model is', x$logLik, 'and the AIC and BIC are', x$AIC,'and', x$BIC, 'respectively.')
  invisible(x)
}
