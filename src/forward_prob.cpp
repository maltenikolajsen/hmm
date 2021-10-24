#include <Rcpp.h>
using namespace Rcpp;

//' forward_prob_cpp
//'
//' @description Calculate matrix of forward probabilities for a given data set and parameters in a HMM
//'
//' @param Gamma A matrix of transition probabilities for the underlying Markov chain
//' @param p A matrix of observation probabilites, i.e. where the i,j'th entry is the likelihood of observing x_j given the hidden state is i
//' @param delta A vector of starting probabilites, i.e. the distribution of the first hidden state
//'
//' @output An m x n matrix of forward probabilities i.e. where the i,j'th entry is the likelihood of observing the first j observations and the j'th hidden state is i
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
        alpha_ji += alpha(k, i-1) * Gamma(k, j) * p(j, i);
        // read as: Prob of being in state k at time i-1 * going from k to j * observing x_i in state j
      }
      alpha(j, i) = alpha_ji;
    }
  }

  return alpha;
}
