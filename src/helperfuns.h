#ifndef HELPERFUNS_H
#define HELPERFUNS_H

#include <RcppArmadillo.h>
#include "tteData.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Function declarations
double logsumexp(arma::vec const& x);
int rcat(arma::vec const& probs);
arma::vec log_stick_breaking(arma::vec const& v);
double rtgamma_cpp(double shape, double rate, double a, double b);
arma::vec stick_breaking(const arma::vec& x, bool log);
arma::vec rdirichlet_rng(int const& n, bool const& log);
arma::mat compute_survival_conditional(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logtimes, bool log);
arma::vec compute_survival_bb(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logomega, arma::vec const& logtimes, bool log);
arma::mat compute_rmst_conditional(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logtimes, bool log);
arma::vec compute_rmst_bb(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logomega, arma::vec const& logtimes, bool log);
arma::mat compute_hazard_conditional(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logtimes, bool log);
arma::vec compute_hazard_bb(arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw, arma::vec const& logomega, arma::vec const& logtimes, bool log);
tteData list2tteData(Rcpp::List const& data);
tteHistData list2tteHistData(Rcpp::List const& histdata);
double rtbeta(double shape1, double shape2, double lower, double upper);
double sample_tbeta_mixture_posterior(double shape1, double shape2, double lower, double upper,double prior_prob, int nexch, int nunexch);
#endif