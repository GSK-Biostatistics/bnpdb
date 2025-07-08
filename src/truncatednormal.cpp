#include "truncatednormal.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

// Exponential rejection sampling (a,inf)
double ers_a_inf(double a) {
  double ainv = 1.0 / a, x, rho;
  do {
    x = R::rexp(ainv) + a;
    rho = std::exp(-0.5 * std::pow(x - a, 2));
  } while (R::runif(0, 1) > rho);
  return x;
}

// Exponential rejection sampling (a,b)
double ers_a_b(double a, double b) {
  double ainv = 1.0 / a, x, rho;
  do {
    x = R::rexp(ainv) + a;
    rho = std::exp(-0.5 * std::pow(x - a, 2));
  } while (R::runif(0, 1) > rho || x > b); 
  return x;
} 

// Normal rejection sampling (a,b)
double nrs_a_b(double a, double b) {
  double x;
  do {
    x = R::rnorm(0, 1);
  } while (x < a || x > b); 
  return x;
}
 
// Normal rejection sampling (a,inf)
double nrs_a_inf(double a) {
  double x;
  do {
    x = R::rnorm(0, 1);
  } while (x < a); 
  return x;
} 

// Half-normal rejection sampling
double hnrs_a_b(double a, double b) {
  double x;
  do {
    x = std::fabs(R::rnorm(0, 1));
  } while (x < a || x > b); 
  return x;
} 

// Uniform rejection sampling
double urs_a_b(double a, double b) {
  double phi_a = R::dnorm(a, 0.0, 1.0, false);
  double ub = (a < 0 && b > 0) ? M_1_SQRT_2PI : phi_a, x;
  do {
    x = R::runif(a, b);
  } while (R::runif(0, 1) * ub > R::dnorm(x, 0, 1, false)); 
  return x;
} 

// Left-truncated normal sampling
double r_lefttruncnorm(double a, double mu, double sigma) {
  double alpha = (a - mu) / sigma;
  return mu + sigma * (alpha < t4 ? nrs_a_inf(alpha) : ers_a_inf(alpha));
} 

// Right-truncated normal sampling
double r_righttruncnorm(double b, double mu, double sigma) {
  return mu - sigma * r_lefttruncnorm(-b, 0.0, 1.0);
} 

// Truncated normal sampling
double r_truncnorm(double a, double b, double mu, double sigma) {
  double alpha = (a - mu) / sigma, beta = (b - mu) / sigma;
  if (beta <= alpha) return NA_REAL;
  double phi_a = R::dnorm(alpha, 0.0, 1.0, false);
  double phi_b = R::dnorm(beta, 0.0, 1.0, false);
  if (alpha <= 0 && 0 <= beta) {
    return mu + sigma * (phi_a <= t1 || phi_b <= t1 ? nrs_a_b(alpha, beta) : urs_a_b(alpha, beta));
  } else if (alpha > 0) { 
    if (phi_a / phi_b <= t2) return mu + sigma * urs_a_b(alpha, beta);
    return mu + sigma * (alpha < t3 ? hnrs_a_b(alpha, beta) : ers_a_b(alpha, beta));
  } else { 
    if (phi_b / phi_a <= t2) return mu - sigma * urs_a_b(-beta, -alpha);
    return mu - sigma * (beta > -t3 ? hnrs_a_b(-beta, -alpha) : ers_a_b(-beta, -alpha));
  }
} 



double rtruncnorm_one_cpp(
  double mu, double sigma, double a, double b
) {
    if (R_FINITE(a) && R_FINITE(b)) {
      return r_truncnorm(a, b, mu, sigma);
    } else if (R_NegInf == a && R_FINITE(b)) {  
      return r_righttruncnorm(b, mu, sigma);
    } else if (R_FINITE(a) && R_PosInf == b) {   
      return r_lefttruncnorm(a, mu, sigma);
    } else if (R_NegInf == a && R_PosInf == b) {   
      return R::rnorm(mu, sigma);
    } 
  return NA_REAL;  // never reached
}

// Vectorized truncated normal sampling
arma::vec rtruncnorm_cpp(
    const arma::vec& mu, const arma::vec& sigma, const arma::vec& a, const arma::vec& b
) {
  int n = a.n_elem;
  arma::vec ret(n);
  for (int i = 0; i < n; ++i) {
    if (R_FINITE(a[i]) && R_FINITE(b[i])) {
      ret[i] = r_truncnorm(a[i], b[i], mu[i], sigma[i]);
    } else if (R_NegInf == a[i] && R_FINITE(b[i])) {
      ret[i] = r_righttruncnorm(b[i], mu[i], sigma[i]);
    } else if (R_FINITE(a[i]) && R_PosInf == b[i]) {
      ret[i] = r_lefttruncnorm(a[i], mu[i], sigma[i]);
    } else if (R_NegInf == a[i] && R_PosInf == b[i]) {
      ret[i] = R::rnorm(mu[i], sigma[i]);
    } else {
      ret[i] = NA_REAL;
    }
  }
  return ret;
}


