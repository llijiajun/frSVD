// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "frpca.h"
#include <Rcpparmadillo.h>
#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
RcppExport SEXP frSVD(arma::sp_mat& X,int k,int q=5) {
  frpca::frPCA fsvd(X,k,q);
  //U=fsvd.matrixU();
  //V=fsvd.matrixV();
  return Rcpp::List::create(
    Rcpp::Named("d")     = fsvd.singularValues(),
    Rcpp::Named("u")     = fsvd.matrixU(),
    Rcpp::Named("v")     = fsvd.matrixV(),
    Rcpp::Named("niter") = (q-1)/2,
    Rcpp::Named("rank")  = k);
}
