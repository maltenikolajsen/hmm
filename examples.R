###################
# EARTHQUAKE DATA #
###################

# Poisson fit
# 2-state
quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
delta <- c(1, 1)/2
lambda <- c(10, 30)

M <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda)

# m-state
m <- 3 # just for example

Gamma <- matrix(0.1, nrow=m, ncol=m) + diag(m) * (1-m * 0.1)
delta <- rep(1/m, m)
lambda <- seq(10, 30, length.out=m)

M <- hmm(quakes, Gamma, delta, dist='poisson', lambda=lambda)

# Binomial fit
quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
delta <- c(1, 1)/2
size <- max(quakes)
prob <- c(10, 30)/size

M <- hmm(quakes, Gamma, delta, dist='binom', size=size, prob=prob)

# Normal fit
quakes <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2
Gamma <- rbind(c(0.9, 0.1), c(0.1, 0.9))
delta <- c(1, 1)/2
mean <- c(10, 30)
sd <- c(1, 1)

M <- hmm(quakes, Gamma, delta, dist='normal', mean=mean, sd=sd)


######################
# GENERATING OWN HMM #
######################

# With a Normal distribution
Gamma <- rbind(c(0.5, 0.25, 0.25),
               c(0.1, 0.8 , 0.1),
               c(  0, 0.2 , 0.8))
delta <- c(1, 0, 0)
mean <- c(0, 5, 10)
sd <- rep(1, 3)

M <- hmm(NULL, Gamma=Gamma, delta=delta, dist='normal', mean=mean, sd=sd)
X <- simulate(M, nsim=100)

# With a custom distribution
Gamma <- rbind(c(0.5, 0.25, 0.25),
               c(0.1, 0.8 , 0.1),
               c(  0, 0.2 , 0.8))
delta <- c(1, 0, 0)
theta <- list(1, 5, 10)
lls <- function(x, param){dunif(x, 0, param)}
lls_mle <- function(x, u){max(x)}
rdist <- function(n, param){runif(n, 0, param)}

M <- hmm(NULL, Gamma=Gamma, delta=delta, lls=lls, param_lls=theta, lls_mle=lls_mle, rdist=rdist)
Z <- simulate(M, nsim=200, include_state=TRUE)
X <- Z[1:200]
G <- Z[201:400]
plot(X, type='h', col=G, lwd=2)

############################
# MIXTURE OF DISTRIBUTIONS #
############################

# Here in state 1 dist is standard normal and in state 2 standard exponential
Gamma <- rbind(c(0.2, 0.8),
               c(0.8, 0.2))
delta <- c(1, 1)/2
param <- list(c(0, 1), 1)

lls <- list(function(x, param){dnorm(x, param[1], param[2])},
            function(x, param){dexp(x, param)})

lls_mle <- list(function(x, u){mean_hat <- sum(u*x) / sum(u); c(mean_hat, sqrt(sum(u*(x-mean_hat)^2) / sum(u)))},
                function(x, u){sum(u)/sum(u*x)})

rdist <- list(function(n, param){do.call(rnorm, args=as.list(c(n, param)))},
              function(n, param){rexp(n, param)})

M <- hmm(NULL, Gamma=Gamma, delta=delta, lls=lls, param_lls=param, lls_mle=lls_mle, rdist=rdist)
Z <- simulate(M, nsim=200, include_state=TRUE)
X <- Z[1:200]
G <- Z[201:400]
plot(X, col=G, pch=19)

# Estimation
M_est <- hmm(X, Gamma=Gamma, delta=delta, lls=lls, param_lls=param, lls_mle=lls_mle, rdist=rdist)
abs(M_est$Gamma-Gamma)
M_est$parameters;param
mean(M_est$viterbi_s == G)
mean(M_est$posterior_s == G)
