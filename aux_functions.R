get_alpha <- function(
  observation, #Observations
  delta, #Starting probabilities
  transition, #Transition matrix
  lls, #Likelihood functions
  param_lls #Parameters for the likelihood functions
){
  m <- length(lls)
  n <- length(observation)
  helper_function <- function(x) {
    emission_prob <- rep(NA, m)
    for(i in 1:m){
      emission_prob[i] <- lls[[i]](x, param_lls[[i]])
    }
    return(emission_prob)
  }

  storage <- as.list(rep(NA, n))

  alpha <- delta * helper_function(observation[1])
  storage[[1]] <- alpha

  lscale <- log(sum(alpha))
  alpha <- alpha / sum(alpha)

  for(i in 2:n){
    alpha <- alpha %*% transition * helper_function(observation[i])
    lscale <- lscale + log(sum(alpha))
    alpha <- alpha / sum(alpha)
    storage[[i]] <- alpha
  }

  storage
}

get_beta <- function(
  observation, #Observations
  delta, #Starting probabilities
  transition, #Transition matrix
  lls, #Likelihood functions
  param_lls #Parameters for the likelihood functions
){

}




#Backwards algorithm


get_beta <- function(
  observaton,
  delta,
  transition,
  lls,
  param_lls
){
  m <- length(lls)
  n <- length(observation)
  helper_function <- function(x) {
    emission_prob <- rep(NA, m)
    for(i in 1:m){
      emission_prob[i] <- lls[[i]](x, param_lls[[i]])
    }
    return(emission_prob)
  }

  beta <- list(rep(NA, n))

  beta[[n]] <- rep(1, m)

  for(i in (n-1):1){
    tmp_beta <- rep(NA, m)
    for(j in seq_along(tmp_beta)){
      tmp_emission_beta <- helper_function(observation[i]) * beta[[i+1]]
      tmp_zeta <- transition[j,] * tmp_emission_beta
      xi_star <- max(tmp_zeta)
      tmp_beta[j] <- xi_star + log(sum(exp(tmp_zeta - xi_star)))
    }
    beta[[i]] <- tmp_beta
  }

  beta
}


get_beta2 <- function(
  observaton,
  delta,
  transition,
  lls,
  param_lls
){
  m <- length(lls)
  n <- length(observation)
  helper_function <- function(x) {
    emission_prob <- rep(NA, m)
    for(i in 1:m){
      emission_prob[i] <- lls[[i]](x, param_lls[[i]])
    }
    return(emission_prob)
  }

  beta <- list(rep(NA, n))

  beta[[n]] <- rep(1, m)

  for(i in (n-1):1){
    tmp_transition_emission <- transition %*% diag(helper_function(observation[i]))
    beta[[i]] <- tmp_transition_emission %*% print(beta[[i+1]])
  }

  beta
}




#EARTHQUAKE DATA
observation <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
m <- 2
lls <- list(
  function(x, param) dpois(x, param[1]),
  function(x, param) dnorm(x, mean = param[1], sd = param[2])#function(x, param) dpois(x, param[1])
)
param_lls <- list(
  c(10),
  c(3, 4)
)
delta <- rep(.5, m)
gamma <- matrix(rep(.5, m^2), nrow = 2)

get_beta(
  observaton,
  delta,
  gamma,
  lls,
  param_lls
)

get_beta2(
  observaton,
  delta,
  gamma,
  lls,
  param_lls
)



#
get_beta3 <- function(
  observaton,
  delta,
  transition,
  lls,
  param_lls
){
  m <- length(lls)
  n <- length(observation)
  helper_function <- function(x) {
    emission_prob <- rep(NA, m)
    for(i in 1:m){
      emission_prob[i] <- lls[[i]](x, param_lls[[i]])
    }
    return(emission_prob)
  }

  logbeta <- list(rep(NA, n))

  logbeta[[n]] <- rep(0, m)

  for(i in (n-1):1){
    tmp_logbeta <- rep(NA, m)
    for(j in seq_along(tmp_logbeta)){
      tmp_zeta <- log(transition[j,]) + log(helper_function(observation[i])) + logbeta[[i+1]]
      logxi_star <- max(tmp_zeta)
      tmp_logbeta[j] <- logxi_star + log(sum(exp(tmp_zeta - logxi_star)))
    }
    logbeta[[i]] <- tmp_logbeta
  }

  beta <- list(rep(NA, n))
  for(i in 1:n)(
    beta[[i]] <- exp(logbeta[[i]])
  )
  beta
}
get_beta3(
  observaton,
  delta,
  gamma,
  lls,
  param_lls
)
#
get_alpha3 <- function(
  observaton,
  delta,
  transition,
  lls,
  param_lls
){
  m <- length(lls)
  n <- length(observation)
  helper_function <- function(x) {
    emission_prob <- rep(NA, m)
    for(i in 1:m){
      emission_prob[i] <- lls[[i]](x, param_lls[[i]])
    }
    return(emission_prob)
  }

  logalpha <- list(rep(NA, n))

  logalpha[[1]] <- log(delta * helper_function(observation[1]))

  for(i in 2:n){
    tmp_logalpha <- rep(NA, m)
    for(j in seq_along(tmp_logalpha)){
      tmp_zeta <- log(transition[j,]) + log(helper_function(observation[i])) + logalpha[[i-1]]
      logxi_star <- max(tmp_zeta)
      tmp_logalpha[j] <- logxi_star + log(sum(exp(tmp_zeta - logxi_star)))
    }
    logalpha[[i]] <- tmp_logalpha
  }

  alpha <- list(rep(NA, n))
  for(i in 1:n)(
    alpha[[i]] <- exp(logalpha[[i]])
  )

  alpha
}
get_alpha3(
  observaton,
  delta,
  gamma,
  lls,
  param_lls
)
