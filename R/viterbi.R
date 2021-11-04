viterbi <- function(obs,delta,gamma,lls,param_lls){
  n <- length(obs); m <- ncol(gamma)

  p <- function(state, x){
    lls[[state]](x, param_lls[[state]])
  }
  p <- Vectorize(p)

  xi <- matrix(NA, ncol = n, nrow = m)
  foo <- delta * p(1:m, obs[1])

  xi[,1] <- foo / sum(foo)
  for(i in 2:n){
    print(xi[,i-1]*gamma)
    foo <- apply(xi[,i-1]*gamma,2,max) * p(1:m, obs[i])
    xi[,i] <- foo/sum(foo)
  }

  iv <- numeric(n)
  iv[n] <- which.max(xi[,n])
  for(i in (n-1):1){
    iv[i] <- which.max(gamma[,iv[i+1]]*xi[,i])
  }
  iv
}

