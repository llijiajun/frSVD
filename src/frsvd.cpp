// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "frpca.h"
#include <Rcpparmadillo.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec frsvd(arma::sp_mat& X,int k,int q=5) {
  frpca::frPCA fsvd(X,k,q);
  //U=fsvd.matrixU();
  //V=fsvd.matrixV();
  arma::vec S=fsvd.singularValues();
  return S;
}
