#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericMatrix forward_prob_cpp(NumericMatrix Gamma, NumericMatrix p) {
  // Gamma = m x m matrix of transition probabilities
  // p = m x n matrix where p(i, j) == P(X_t=x_j | C_t=i)
  
  // Set n and m for sake of notation
  int n = p.ncol();
  int m = p.nrow();
  
  // Initialize matrix of backward probabilities
  NumericMatrix beta(m, n);
  
  // Set last column of beta
  for(int i=0; i<m; i++){
    beta(i, n-1) = 1;
  }
  
  // Recursively fill out the rest of beta
  for(int i=n-1; i>0; i--){
    for(int j=0; j<m; j++){
      double beta_ji = 0;
      for(int k=0; k<m; k++){
        beta_ji += Gamma(j, k) * p(k, i) * beta(j, i);
      }
      beta(j, i-1) = beta_ji;
    }
  }
  
  return beta;
}