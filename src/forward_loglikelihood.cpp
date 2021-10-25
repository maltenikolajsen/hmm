#include <Rcpp.h>
using namespace Rcpp;

//' forward_ll_cpp
//'
//' @description Calculate matrix of forward log-probabilities for a given data set and parameters in a HMM
//'
//' @param Gamma A matrix of transition log-probabilities for the underlying Markov chain
//' @param p A matrix of observation log-probabilites, i.e. where the i,j'th entry is the log-likelihood of observing x_j given the hidden state is i
//' @param delta A vector of starting log-probabilites, i.e. the log of the distribution of the first hidden state
//'
//' @output An m x n matrix of forward log-probabilities i.e. where the i,j'th entry is the likelihood of observing the first j observations and the j'th hidden state is i
// [[Rcpp::export]]
NumericMatrix forward_ll_cpp(NumericMatrix Gamma, NumericMatrix p, NumericVector delta) {
  // Set n and m for sake of notation
  int n = p.ncol();
  int m = p.nrow();

  // Initialize matrix of forward log-probabilities
  NumericMatrix logalpha(m, n);

  // Set first column of log-alpha
  for(int i=0; i<m; i++){
    logalpha(i, 0) = delta[i] + p(i, 0);
  }

  // Recursively fill out the rest of log-alpha
  for(int i=1; i<n; i++){
    double c = max(logalpha(_, i-1));
    for(int j=0; j<m; j++){
      double logalpha_ji = 0;
      for(int k=0; k<m; k++){
        logalpha_ji += exp(logalpha(k, i-1) - c + Gamma(k, j) + p(j, i));
      }
      logalpha(j, i) = c+log(logalpha_ji);
    }
  }

  return logalpha;
}
