fa.ll <- function(
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

  alpha <- delta * helper_function(observation[1])
  lscale <- log(sum(alpha))
  alpha <- alpha / sum(alpha)

  for(i in 2:n){
    alpha <- alpha %*% transition * helper_function(observation[i])
    lscale <- lscale + log(sum(alpha))
    alpha <- alpha / sum(alpha)
  }

  -lscale
}

hmm.pn2pw <- function(m, delta, gamma)
{
  foo     <- log(gamma/diag(gamma))
  tgamma  <- as.vector(foo[!diag(m)])
  tdelta <- log(delta[-1]/delta[1])
  parvect <- c(tdelta,tgamma)
  return(parvect)
}

hmm.pw2pn <- function(m, parvect){
  {
    k <- length(parvect)

    foo <- c(1, exp(parvect[1:(m-1)]))
    delta <- foo/sum(foo)

    gamma         <- diag(m)
    gamma[!gamma] <- exp(parvect[m:k])
    gamma         <- gamma/apply(gamma,1,sum)

    delta<-foo/sum(foo)}
  return(list(gamma=gamma,delta=delta))
}


hmm.fac.obj_param_to_neutral_param <- function(m, param_lls) {
  hmm.obj_param_to_neutral_param <- function(
    param #Whack ass parameters
  )
  {
    delta_index <- 1:(m-1)
    gamma_index <- m:(m^2-1)
    tmp <- hmm.pw2pn(m, param[1:(m^2-1)])

    start_index_param <- m^2
    remaining_param <- param[start_index_param:length(param)]

    tmp_param <- param_lls
    k <- 1
    while(k <= length(remaining_param)){
      for(i1 in seq_along(tmp_param)){
        for(i2 in seq_along(tmp_param[[i1]])){
          tmp_param[[i1]][i2] <- remaining_param[k]
          k <- k+1
        }
      }
    }
    return(
      list(
        "delta" = tmp$delta,
        "gamma" = tmp$gamma,
        "param_lls" = tmp_param
      )
    )
  }
}

fa.ll.wrapper <- function(observation, m, lls, param_lls){
  obj <- function(param, w2n){
    input_list <- w2n(param)
    print(input_list)
    fa.ll(observation,
          input_list$delta,
          input_list$gamma,
          lls,
          input_list$param_lls)
  }
  obj
}


hehexd <- function(
  observation,
  m,
  lls,
  param_lls,
  param_lls_lb,
  param_lls_ub,
  delta,
  gamma,
  no_iter = 1000
){
  no_param <- length(param_lls_ub)
  w2n <- hmm.fac.obj_param_to_neutral_param(m, param_lls)
  my_obj <- fa.ll.wrapper(observation, m, lls)
  wp <- c(hmm.pn2pw(m, delta, gamma), as.vector(sapply(param_lls, function(x) x)))
  obj_foo <- function(x) my_obj(x, w2n)
  obj_foo(wp)
  lb <- c(rep(-Inf, m-1), rep(-Inf, m^2-m), param_lls_lb)
  ub <- c(rep(Inf, m-1), rep(Inf, m^2-m), param_lls_ub)

  res1 <- nloptr::nloptr(
    x0 = wp,
    eval_f = obj_foo,
    lb = lb,
    ub = ub,
    opts = list("algorithm" ="NLOPT_LN_COBYLA",
                "xtol_rel" = 1.0e-15,
                "maxeval" = no_iter)
  )
  w2n(res1$solution)
}

#Stupid data
observation <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
m <- 2
lls <- list(
  function(x, param) dnorm(x, mean = param[1], sd = param[2]),
  function(x, param) dnorm(x, mean = param[1], sd = param[2])
)
param_lls <- list(
  c(1, 2),
  c(3, 4)
)
param_lls_lb <- c(-Inf, 0, -Inf, 0)
param_lls_ub <- c(Inf, Inf, Inf, Inf)
delta <- rep(.5, m)
gamma <- matrix(rep(.5, m^2), nrow = 2)

hehexd(observation,
       m,
       lls,
       param_lls,
       param_lls_lb,
       param_lls_ub,
       delta,
       gamma,
       no_iter = 5000)

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
param_lls_lb <- c(0,-Inf,0)
param_lls_ub <- c(Inf, Inf, Inf)
delta <- rep(.5, m)
gamma <- matrix(rep(.5, m^2), nrow = 2)
a <- hehexd(observation,
            m,
            lls,
            param_lls,
            param_lls_lb,
            param_lls_ub,
            delta,
            gamma,
            no_iter = 5000)
a


m <- 3
lls <- list(
  function(x, param) dpois(x, param[1]),
  function(x, param) dpois(x, param[1]),
  function(x, param) dpois(x, param[1])
)
param_lls <- list(
  c(10),
  c(10),
  c(100)
)
param_lls_lb <- c(0,0,0)
param_lls_ub <- c(Inf, Inf, Inf)
delta <- rep(1/m, m)
gamma <- matrix(rep(1/m, m^2), nrow = 2)
b <- hehexd(observation,
            m,
            lls,
            param_lls,
            param_lls_lb,
            param_lls_ub,
            delta,
            gamma,
            no_iter = 5000)
b
