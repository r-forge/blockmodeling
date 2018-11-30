#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::NumericVector test( const arma::Cube<double> M, const Rcpp::NumericVector clu, const arma::Cube<double> weights, const Rcpp::NumericVector n, const Rcpp::NumericVector nClu )
{
    return clu;
}

