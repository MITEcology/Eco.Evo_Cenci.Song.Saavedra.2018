#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <iostream>     // std::cout
#include <valarray>     // std::valarray
#include <Rcpp.h>

// Nota bene, Armadillo does not change the matrix then there is no need to clone
// Otherwise then remember to clone


// [[Rcpp::export]]
arma::vec FixedPoint(arma::mat X, arma::vec Y){
  /// This function return the abundances of a linear system
  return(inv(X) * Y);
}

// [[Rcpp::export]]
arma::vec getEigenvalues(arma::mat M) {
  /// Return the real part of the eigenvalues
  return real(eig_gen(M));
}

// [[Rcpp::export]]
arma::uword check_if_stable(arma::mat M) {
  /// This function return:
  /// 0 if there are negative eigenvalues
  /// 1 otherwise
  arma::vec autovalori = real(eig_gen(M));
  arma::uvec ids = find(autovalori <= 0); // Find indices
  if(ids.size() > 0){
    return 0;
  }else{
    return(1);
  }
}

// [[Rcpp::export]]
arma::colvec Centroid(arma::mat X){
  return sum(X,1)/X.n_cols;
}
