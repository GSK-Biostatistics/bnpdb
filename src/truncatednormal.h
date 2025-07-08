
#ifndef TRUNCATEDNORMAL_H
#define TRUNCATEDNORMAL_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double rtruncnorm_one_cpp(double mu, double sigma, double a, double b);
arma::vec rtruncnorm_cpp(const arma::vec& mu, const arma::vec& sigma, const arma::vec& a, const arma::vec& b);
#endif
