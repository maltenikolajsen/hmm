#include <Rcpp.h>
using namespace Rcpp;

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
    if(c == R_NegInf){
      c = 0;
    }
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
