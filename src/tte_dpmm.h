#ifndef TTE_DPMM_H
#define TTE_DPMM_H

#include <RcppArmadillo.h>
#include "truncatednormal.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace TTE_DPMM {

// BaseMeasure class definition
  class BaseMeasure {
  public: 
    arma::vec beta_mean;
    arma::mat beta_prec;
    double tau_shape;
    double tau_rate;
    arma::mat beta_cov;
    arma::vec beta_prec_mean;
    double beta_mean_prec_mean;
    double alpha_shape;
    double alpha_rate;
    double alpha_lower;
    double alpha_upper;
    
    BaseMeasure(
      const arma::vec& beta_mean_, const arma::mat& beta_prec_,
      double tau_shape_, double tau_rate_,
      double alpha_shape_, double alpha_rate_, double alpha_lower_, double alpha_upper_
    ) {
      // Exchangeable
      beta_mean      = beta_mean_;
      beta_prec      = beta_prec_;
      tau_shape      = tau_shape_;
      tau_rate       = tau_rate_;
      alpha_shape    = alpha_shape_;
      alpha_rate     = alpha_rate_;
      alpha_lower    = alpha_lower_;
      alpha_upper    = alpha_upper_;
      beta_cov       = arma::inv_sympd(beta_prec);
      beta_prec_mean      = beta_prec * beta_mean;
      beta_mean_prec_mean = arma::dot(beta_mean, beta_prec * beta_mean);
    }  
  };   // BaseMeasure
  
  // Pars class definition
  class Pars {
  public: 
    // Exchangeable Pars
    int K;                  // Number of classes (exchangeable)
    double alpha;           // DP concentration parameter
    arma::mat beta;         // List of beta vectors (one per class)
    arma::vec tau;          // Precision values (shared across classes)
    arma::vec sigma;        // standard deviation (as opposed to precision)
    arma::vec v;            // Stick-breaking weights
    arma::vec logw;         // Class probabilities
    arma::uvec z;           // Latent class indicators
    arma::vec ystar;        // Latent variable for time-to-event data
    arma::uvec sizes;       // Sizes of each class
    arma::mat Mu;           // Matrix of means for current data
    
    // Constructor
    Pars(
      tteData const& data
    // Exchangeable
    , arma::mat const& beta_init, arma::vec const& tau_init
    , double alpha_init, arma::vec const& v_init
    , arma::uvec const& z_init
    ) {
      // Exchangeable
      alpha = alpha_init;
      beta  = beta_init;
      K     = beta.n_cols;
      tau   = tau_init;
      sigma = 1.0 / arma::sqrt(tau);
      v     = v_init;
      logw  = stick_breaking(v, true);                     // Compute log probabilities
      z     = z_init;
      Mu    = data.X * beta;          // Initialize means
      ystar = data.y;
      ystar.elem(data.cenindx) += 0.1; // Adjust for censored
      sizes = arma::zeros<arma::uvec>(K);
      for ( int i = 0; i < data.n; i++ )
        sizes[z[i]] += 1;
    }
    
    // Update regression parameters
    void update_regression_params(const tteData& data, const BaseMeasure& bm) {
      for (int k = 0; k < K; k++) {
        arma::uvec class_indices_cur  = arma::find(z == k); // indices of current data subjects
        int n_k  = class_indices_cur.n_elem;                // total number of subjects in class
        if (n_k == 0) {
          tau[k]      = R::rgamma(bm.tau_shape, 1.0 / bm.tau_rate);
          beta.col(k) = arma::mvnrnd(bm.beta_mean, bm.beta_cov);
        } else {
          arma::mat X_k       = data.X.rows(class_indices_cur);  // class-specific design matrix
          arma::vec ystar_k   = ystar(class_indices_cur);        // class-specific outcome
          arma::mat XtX       = X_k.t() * X_k;
          arma::vec Xty       = X_k.t() * ystar_k;
          
          // Update beta | tau
          arma::mat beta_prec = bm.beta_prec + tau(k) * XtX;
          arma::mat beta_cov  = arma::inv_sympd(beta_prec);
          arma::vec beta_mean = beta_cov * (tau(k) * Xty + bm.beta_prec_mean);
          beta.col(k)         = arma::mvnrnd(beta_mean, beta_cov);
          
          // Update tau | beta
          double sse       = arma::sum( arma::square( ystar_k - X_k * beta.col(k) ) );
          double tau_shape = bm.tau_shape + 0.5 * n_k;
          double tau_rate  = bm.tau_rate + 0.5 * sse;
          tau[k]           = R::rgamma(tau_shape, 1.0 / tau_rate);
        }
      }
      // Update means and standard deviation
      sigma    = 1.0 / arma::sqrt(tau);
      Mu       = data.X * beta;
    }  // update_regression_params
    
    // Update class assignments (current)
    void update_classes(const tteData& data) {
      int n        = data.n;
      arma::vec logprobs(K, arma::fill::zeros);
      arma::vec probs(K, arma::fill::zeros);
      for (int i = 0; i < n; i++) {
        sizes[z[i]] -= 1;
        logprobs = logw + arma::log_normpdf(ystar[i], Mu.row(i).t(), sigma);
        probs    = arma::exp(logprobs - logsumexp(logprobs));
        z[i]     = rcat(probs);
        sizes[z[i]] += 1;
      }
    }   // update_classes
    
    // Update stick-breaking weights (exchangeable)
    void update_sb_variables() {
      double higher_sizes = sizes[K - 1];
      for (int k = K - 2; k >= 0; k--) {
        double post_shape1 = sizes[k] + 1;
        double post_shape = alpha + higher_sizes;
        v[k]          = R::rbeta(post_shape1, post_shape);
        higher_sizes += sizes[k];
      }
      logw = stick_breaking(v, true);
    }  // update_sb_variables
    
    // Update concentration parameter (exchangeable)
    void update_concentration_parameter(const BaseMeasure& bm) {
      double post_shape = bm.alpha_shape + K - 1;
      double post_rate  = bm.alpha_rate - arma::sum(arma::log1p(-v));
      alpha = rtgamma_cpp(post_shape, post_rate, bm.alpha_lower, bm.alpha_upper);
    }  // update_concentration_parameter
    
    // Update censored variables (current)
    void update_censored_vars(const tteData& data) {
      if (data.ncen > 0) {
        for (int j = 0; j < data.ncen; j++) {
          int i   = data.cenindx[j];   // index of censored value
          int z_i = z[i];              // class of censored value
          double mean_i    = Mu(i, z_i);
          double sigma_i   = sigma(z_i);
          ystar[i] = rtruncnorm_one_cpp(mean_i, sigma_i, data.y(i), R_PosInf);
        }
      }
    }  // update_censored_vars
  }; // Pars
  
  
  
  // Convert Rcpp list of parameters to Pars object
  Pars list2pars(
      tteData const& data, Rcpp::List const& pars
  ) {
    arma::mat beta  = pars["beta"];
    arma::vec tau   = pars["tau"];
    double alpha    = pars["alpha"];
    arma::vec v     = pars["v"];
    arma::uvec z    = pars["z"];
    Pars params(data, beta, tau, alpha, v, z);
    return params;
  }
  
  // Convert Rcpp list of base measure values to basemeasure object
  BaseMeasure list2basemeasure(Rcpp::List const& basemeasure) {
    // Base measure for exchangeable parameters
    arma::vec beta_mean = basemeasure["beta_mean"];
    arma::mat beta_prec = basemeasure["beta_prec"];
    double tau_shape  = basemeasure["tau_shape"];
    double tau_rate = basemeasure["tau_rate"];
    double alpha_shape = basemeasure["alpha_shape"];
    double alpha_rate = basemeasure["alpha_rate"];
    double alpha_lower = basemeasure["alpha_lower"];
    double alpha_upper = basemeasure["alpha_upper"];
    BaseMeasure bm(
        beta_mean, beta_prec, tau_shape, tau_rate
      , alpha_shape, alpha_rate, alpha_lower, alpha_upper 
    );
    return bm;
  }
} // namespace TTE_DPMM

#endif // TTE_DPMM_H
