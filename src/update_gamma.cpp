#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_log_gamma_cpp(NumericMatrix log_alpha, NumericMatrix log_beta, NumericMatrix log_gamma, NumericMatrix log_p_mat, float log_ll) {
  // Initialize n and m for notation
  int n = log_alpha.ncol();
  int m = log_alpha.nrow();

  // First calculate log(n_{i, j}) by log-sum-exp trick
  NumericMatrix log_n(m, m);
  for(int i=0; i<m; i++){
    for(int j=0; i<m; j++){
      NumericVector log_v_ij(n);
      for(int k=1; k<n; k++){
        log_v_ij[k] = log_alpha(i, k-1) + log_gamma(i, j) + log_p_mat(j, k) + log_beta(j, k) - log_ll;
      }
      double C_1 = max(log_v_ij[Range(1, n-1)]);
      if(C_1 == R_NegInf){
        C_1 = 0;
      }
      double n_ij = 0;
      for(int k=1; k<n; k++){
        n_ij += exp(log_v_ij(k) - C_1);
      }
      log_n(i, j) = C_1 + log(n_ij);
    }
  }

  // Next calculate log(Gamma_{i, j}) by log-sum-exp trick
  NumericMatrix log_Gamma(m, m);
  for(int i=0; i<m; i++){
    for(int j=0; j<m; j++){
      double C_2 = max(log_n(i, _ ));
      double n_ik = 0;
      for(int k=0; k<m; k++){
        n_ik += exp(log_n(i, k) - C_2);
      }
      log_Gamma(i, j) = log_n(i, j) - (C_2 + log(n_ik));
    }
  }

  return log_Gamma;
}
