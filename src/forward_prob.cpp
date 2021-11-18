#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix forward_prob_cpp(NumericMatrix Gamma, NumericMatrix p, NumericVector delta) {
  // Set n and m for sake of notation
  int n = p.ncol();
  int m = delta.size();

  // Initialize matrix of forward probabilities
  NumericMatrix alpha(m, n);

  // Set first column of alpha
  for(int i=0; i<m; i++){
    alpha(i, 0) = delta[i] * p(i, 0);
  }

  // Recursively fill out the rest of alpha
  for(int i=1; i<n; i++){
    for(int j=0; j<m; j++){
      double alpha_ji = 0;
      for(int k=0; k<m; k++){
        alpha_ji += alpha(k, i-1) * Gamma(k, j);
      }
      alpha(j, i) = alpha_ji * p(j, i);
    }
  }

  return alpha;
}
