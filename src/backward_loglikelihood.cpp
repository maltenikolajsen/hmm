#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix backward_ll_cpp(NumericMatrix Gamma, NumericMatrix p) {
  // Set n and m for sake of notation
  int n = p.ncol();
  int m = p.nrow();

  // Initialize matrix of backward log-probabilities
  NumericMatrix logbeta(m, n);

  // Set last column of log-beta
  for(int i=0; i<m; i++){
    logbeta(i, n-1) = 0;
  }

  // Recursively fill out the rest of beta
  for(int i=n-1; i>0; i--){
    double c = max(logbeta(_, i));
    if(c == R_NegInf){
      c = 0;
    }
    for(int j=0; j<m; j++){
      double beta_ji = 0;
      for(int k=0; k<m; k++){
        beta_ji += exp(Gamma(j, k) + p(k, i) + logbeta(k, i) - c);
      }
      logbeta(j, i-1) = c + log(beta_ji);
    }
  }

  return logbeta;
}
