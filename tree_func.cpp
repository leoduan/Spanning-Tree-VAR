#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int gumbelMaxC(const vec& logWeight){
  int n = logWeight.n_elem;
  vec noisy_weights = logWeight + (-log(-log(randu(n))));
  int v = noisy_weights.index_max();
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat randomWalkCover(const mat& logW) {
  // Initialize variables
  int n = logW.n_rows;
  arma::mat A_T_(n, n, arma::fill::zeros);
  arma::Col<int> InTree(n, arma::fill::zeros);
  arma::Col<int> Next(n);
  
  // Set up Next and InTree
  int r = 0;
  InTree(r) = 1;
  Next(r) = 0;
  
  // Compute A_T_
  for (int i = 0; i < n; i++) {
    int u = i;
    //do a random walk, until getting back to the tree
    while (!InTree(u)) {
      int v = gumbelMaxC(logW.col(u));
      Next(u) = v;
      u = v;
    }
    u = i;
    //go through this path that we just walked, change the label of each node to InTree = True
    while (!InTree(u)) {
      InTree(u) = 1;
      u = Next(u);
    }
  }
  
  // Construct adjacency matrix
  for (int u = 1; u < n; u++) {
    A_T_(u, Next(u)) = 1;
    A_T_(Next(u), u) = 1;
  }
  
  return A_T_;
}



