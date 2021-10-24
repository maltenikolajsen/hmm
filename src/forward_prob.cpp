#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericMatrix forward_prob_cpp(NumericMatrix Gamma, NumericMatrix p, NumericVector delta) {
  // Gamma = m x m matrix of transition probabilities
  // p = m x n matrix where p(i, j) == P(X_t=x_j | C_t=i)
  // delta = vector of length m of starting probabilities
  
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
        alpha_ji += alpha(k, i-1) * Gamma(k, j) * p(j, i); // read as: Prob of being in state k at time i-1 * going from k to j * observing x_i in state j
      }
      alpha(j, i) = alpha_ji;
    }
  }
  
  return alpha;
}