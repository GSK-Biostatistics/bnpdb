#ifndef TTE_BNPDB_H
#define TTE_BNPDB_H

#include <RcppArmadillo.h>
#include "truncatednormal.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace TTE_BNPDB {
  // BaseMeasure class definition
  class BaseMeasure {
  public: 
    // Exchangeable parameters
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
    // Unexchangeable parameters
    arma::vec beta_unexch_mean;
    arma::mat beta_unexch_prec;
    double tau_unexch_shape;
    double tau_unexch_rate;
    arma::mat beta_unexch_cov;
    arma::vec beta_unexch_prec_mean;
    double beta_unexch_mean_prec_mean;
    double alpha_unexch_shape;
    double alpha_unexch_rate;
    double alpha_unexch_lower;
    double alpha_unexch_upper;
    double pexch_shape1;
    double pexch_shape2;
    double pexch_lower;
    double pexch_upper;
    double pexch_priorprob;
    
    BaseMeasure(
      const arma::vec& beta_mean_, const arma::mat& beta_prec_,
      double tau_shape_, double tau_rate_,
      double alpha_shape_, double alpha_rate_, double alpha_lower_, double alpha_upper_
    , const arma::vec& beta_unexch_mean_, const arma::mat& beta_unexch_prec_
    , double tau_unexch_shape_, double tau_unexch_rate_
    , double alpha_unexch_shape_, double alpha_unexch_rate_, double alpha_unexch_lower_, double alpha_unexch_upper_
    , double pexch_shape1_, double pexch_shape2_, double pexch_lower_, double pexch_upper_, double pexch_priorprob_
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
      
      // Unexchangeable
      beta_unexch_mean   = beta_unexch_mean_;
      beta_unexch_prec   = beta_unexch_prec_;
      tau_unexch_shape   = tau_unexch_shape_;
      tau_unexch_rate    = tau_unexch_rate_;
      alpha_unexch_shape = alpha_unexch_shape_;
      alpha_unexch_rate  = alpha_unexch_rate_;
      alpha_unexch_lower = alpha_unexch_lower_;
      alpha_unexch_upper = alpha_unexch_upper_;
      beta_unexch_cov    = arma::inv_sympd(beta_unexch_prec);
      beta_unexch_prec_mean      = beta_unexch_prec * beta_unexch_mean;
      beta_unexch_mean_prec_mean = arma::dot(beta_unexch_mean, beta_unexch_prec * beta_unexch_mean);
      pexch_shape1         = pexch_shape1_;
      pexch_shape2         = pexch_shape2_;
      pexch_lower          = pexch_lower_;
      pexch_upper          = pexch_upper_;
      pexch_priorprob      = pexch_priorprob_;
    }  
  };
  
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
    
    // Unxchangeable Pars
    int K0;                  // Number of classes (unexchangeable)
    double alpha_unexch;           // DP concentration parameter
    arma::mat beta_unexch;         // List of beta_unexch vectors (one per class)
    arma::vec tau_unexch;          // Precision values (shared across classes)
    arma::vec sigma0;        // standard deviation (as opposed to precision)
    arma::vec v0;            // Stick-breaking weights
    arma::vec logw0;         // Class probabilities
    arma::uvec z0;           // Latent class indicators
    arma::vec ystar0;        // Latent variable for time-to-event data
    arma::uvec sizes0;       // sizes0 of each class
    arma::mat Mu0;           // Matrix of means for current data
    double pexch;            // probability of being exchangeable
    arma::uvec eps0;         // 1 = exchangeable; 0 = unexchangeable
    arma::mat Mu0_exch;      // Matrix of exchangeable means
    arma::mat Mu0_unexch;    // Matrix of unexchangeable means
    int n0exch;              // number of exchangeable subjects
    int n0unexch;            // number of unexchangeable subjects
    
    // Constructor
    Pars(
      tteData const& data, tteHistData const& histdata
    // Exchangeable
    , arma::mat const& beta_init, arma::vec const& tau_init
    , double alpha_init, arma::vec const& v_init
    , arma::uvec const& z_init
    // Unexchangeable
    , arma::mat const& beta_unexch_init, arma::vec const& tau_unexch_init
    , double alpha_unexch_init, arma::vec const& v0_init
    , double const& pexch_init
    , arma::uvec const& eps0_init, arma::uvec const& z0_init
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
      
      // Unexchangeable
      alpha_unexch = alpha_unexch_init;
      beta_unexch  = beta_unexch_init;
      K0     = beta_unexch.n_cols;
      tau_unexch   = tau_unexch_init;
      sigma0 = 1.0 / arma::sqrt(tau_unexch);
      v0     = v0_init;
      logw0  = stick_breaking(v0, true);                     // Compute log probabilities
      z0     = z0_init;
      ystar0 = histdata.y;
      ystar0.elem(histdata.cenindx) += 0.1; // Adjust for censored
      pexch      = pexch_init;
      Mu0_exch   = histdata.X * beta;
      Mu0_unexch = histdata.Xunexch * beta_unexch;
      eps0  = eps0_init;
      sizes0 = arma::zeros<arma::uvec>(K0);
      n0exch   = 0;
      n0unexch = 0;
      for ( int i = 0; i < histdata.n; i++ ) {
        if ( eps0[i] == 1 ) {
          sizes[z0[i]] += 1;
          n0exch += 1;
        }
        else {
          sizes0[z0[i]] += 1;
          n0unexch      += 1;
        }
      }
    }
    
    // Update regression parameters (exchangeable)
    void update_regression_params_exch(const tteData& data, const tteHistData& histdata, const BaseMeasure& bm) {
      for (int k = 0; k < K; k++) {
        arma::uvec class_indices_cur  = arma::find(z == k);               // indices of current data subjects
        arma::uvec class_indices_hist = arma::find(z0 == k && eps0 == 1);  // indices of historical data subjects
        int N_k  = class_indices_cur.n_elem;         // number of current data subjects in class k
        int N0_k = class_indices_hist.n_elem;        // number of current data subjects in class k
        int n_k  = N_k + N0_k;                       // total number of subjects in class (current + historical)
        if (n_k == 0) {
          tau[k]      = R::rgamma(bm.tau_shape, 1.0 / bm.tau_rate);
          beta.col(k) = arma::mvnrnd(bm.beta_mean, bm.beta_cov);
        } else {
          arma::mat X_k(n_k, data.p, arma::fill::zeros);
          arma::vec ystar_k(n_k, arma::fill::zeros);
          // Fill first N_k elements with current data assigned to class k
          if ( N_k > 0 ) {
            X_k.head_rows(N_k)  = data.X.rows(class_indices_cur);
            ystar_k.head(N_k)   = ystar(class_indices_cur);
          }
          // Fill remaining N0_k elmements with exchangeable historical data
          // assigned to class k
          if ( N0_k > 0 ) {
            X_k.tail_rows(N0_k) = histdata.X.rows(class_indices_hist);
            ystar_k.tail(N0_k)  = ystar0(class_indices_hist);
          }
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
      Mu0_exch = histdata.X * beta;
    }
    
    
    void update_regression_params_unexch(const tteHistData& histdata, const BaseMeasure& bm) {
      for (int k = 0; k < K0; k++) {
        arma::uvec class_indices_hist = arma::find(z0 == k && eps0 == 0);  // indices of historical data subjects
        int n_k = class_indices_hist.n_elem;        // total number of subjects in class
        if (n_k == 0) {
          tau_unexch[k]      = R::rgamma(bm.tau_unexch_shape, 1.0 / bm.tau_unexch_rate);
          beta_unexch.col(k) = arma::mvnrnd(bm.beta_unexch_mean, bm.beta_unexch_cov);
        } else { 
          arma::mat X_k     = histdata.Xunexch.rows(class_indices_hist);
          arma::vec ystar_k = ystar0(class_indices_hist);
          arma::mat XtX     = X_k.t() * X_k;
          arma::vec Xty     = X_k.t() * ystar_k;
          
          // Update beta_unexch | tau_unexch
          arma::mat beta_unexch_prec = bm.beta_unexch_prec + tau_unexch(k) * XtX;
          arma::mat beta_unexch_cov  = arma::inv_sympd(beta_unexch_prec);
          arma::vec beta_unexch_mean = beta_unexch_cov * (tau_unexch(k) * Xty + bm.beta_unexch_prec_mean);
          beta_unexch.col(k)         = arma::mvnrnd(beta_unexch_mean, beta_unexch_cov);
          
          // Update tau_unexch | beta_unexch
          double sse          = arma::sum( arma::square( ystar_k - X_k * beta_unexch.col(k) ) );
          double tau_shape    = bm.tau_shape + 0.5 * n_k;
          double tau_rate     = bm.tau_rate  + 0.5 * sse;
          tau_unexch[k]       = R::rgamma(tau_shape, 1.0 / tau_rate);
        } 
      } 
      // Update means and standard deviation
      sigma0     = 1.0 / arma::sqrt(tau_unexch);
      Mu0_unexch = histdata.Xunexch * beta_unexch;
    } 
    
    // Update class assignments (current)
    void update_classes_current(const tteData& data) {
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
    }
    
    void update_classes_historical(const tteHistData& histdata) {
      int n0 = histdata.n;
      arma::vec logprobs(K + K0, arma::fill::zeros);
      arma::vec probs(K + K0, arma::fill::zeros);
      double logp   = log(pexch);
      double log1mp = log1p(-pexch);
      int z0temp = 0;
      int z0_i = 0;
      int eps0_i = 0;
      for ( int i = 0; i < n0; i++ ) {
        z0_i   = z0[i];
        eps0_i = eps0[i];
        if ( eps0_i ) {
          sizes[z0_i] -= 1;
          n0exch      -= 1;
        }
        else {
          sizes0[z0_i] -= 1;
          n0unexch     -= 1;
        }
        
        logprobs.head(K)   = logp   + logw  + arma::log_normpdf(ystar0[i], Mu0_exch.row(i).t(), sigma);
        logprobs.tail(K0)  = log1mp + logw0 + arma::log_normpdf(ystar0[i], Mu0_unexch.row(i).t(), sigma0);
        logprobs          -= logprobs.max();
        logprobs          -= std::log(arma::sum(arma::exp(logprobs)));
        probs              = exp(logprobs);
        z0temp             = rcat(probs);
        
        // Get new class
        if ( z0temp < K ) {
          // Subject is exchangeable
          z0[i]    = z0temp;
          eps0[i]  = 1;
          n0exch  += 1;
          sizes[z0[i]] += 1;
        }
        else {
          z0[i]     = z0temp - K;
          eps0[i]   = 0;
          n0unexch += 1;
          sizes0[z0[i]] += 1;
        }
      }
    } 
    
    // Update stick-breaking weights (exchangeable)
    void update_sb_variables_exch() {
      double higher_sizes = sizes[K - 1];
      for (int k = K - 2; k >= 0; k--) {
        double post_shape1 = sizes[k] + 1;
        double post_shape2 = alpha + higher_sizes;
        v[k]          = R::rbeta(post_shape1, post_shape2);
        higher_sizes += sizes[k];
      }
      logw = stick_breaking(v, true);
    }
    
    // Update stick-breaking weights (unexchangeable)
    void update_sb_variables_unexch() {
      double higher_sizes = sizes0[K0 - 1];
      for (int k = K0 - 2; k >= 0; k--) {
        double post_shape1 = sizes0[k] + 1;
        double post_shape2 = alpha_unexch + higher_sizes;
        v0[k] = R::rbeta(post_shape1, post_shape2);
        higher_sizes += sizes0[k];
      }
      logw0 = stick_breaking(v0, true);
    }
    
    // Update concentration parameter (exchangeable)
    void update_concentration_parameter_exch(const BaseMeasure& bm) {
      double post_shape = bm.alpha_shape + K - 1;
      double post_rate  = bm.alpha_rate - arma::sum(arma::log1p(-v));
      alpha = rtgamma_cpp(post_shape, post_rate, bm.alpha_lower, bm.alpha_upper);
    }
    
    // Update concentration parameter (unexchangeable)
    void update_concentration_parameter_unexch(const BaseMeasure& bm) {
      double post_shape = bm.alpha_unexch_shape + K0 - 1;
      double post_rate  = bm.alpha_unexch_rate - arma::sum(arma::log1p(-v0));
      alpha_unexch  = rtgamma_cpp(post_shape, post_rate, bm.alpha_unexch_lower, bm.alpha_unexch_upper);
    }
    
    // Update censored variables (current)
    void update_censored_vars_current(const tteData& data) {
      if (data.ncen > 0) {
        for (int j = 0; j < data.ncen; j++) {
          int i   = data.cenindx[j];   // index of censored value
          int z_i = z[i];              // class of censored value
          double mean_i    = Mu(i, z_i);
          double sigma_i   = sigma(z_i);
          ystar[i] = rtruncnorm_one_cpp(mean_i, sigma_i, data.y(i), R_PosInf);
        }
      }
    }  
    
    // Update censored variables (historical)
    void update_censored_vars_historical(const tteHistData& histdata) {
      if (histdata.ncen > 0) {
        for (int j = 0; j < histdata.ncen; j++) {
          int i      = histdata.cenindx[j];   // index of censored value
          int eps0_i = eps0[i];               // exchangeability indicator
          int z0_i   = z0[i];                 // class of censored value
          double mean_i;
          double sigma_i;
          if ( eps0_i ) {
            mean_i  = Mu0_exch(i, z0_i);
            sigma_i = sigma[z0_i];
          } 
          else {
            mean_i  = Mu0_unexch(i, z0_i);
            sigma_i = sigma0[z0_i];
          }
          ystar0[i] = rtruncnorm_one_cpp(mean_i, sigma_i, histdata.y(i), R_PosInf);
        }
      }
    }
    
    void update_pexch(BaseMeasure const& bm) {
      pexch = sample_tbeta_mixture_posterior(
        bm.pexch_shape1, bm.pexch_shape2, bm.pexch_lower, bm.pexch_upper
      , bm.pexch_priorprob, n0exch, n0unexch
      );
    }
  };



  // Convert Rcpp list of parameters to Pars object
  Pars list2pars(
      tteData const& data, tteHistData const& histdata, Rcpp::List const& pars
  ) {
    // Exchangeable params
    arma::mat beta  = pars["beta"];
    arma::vec tau   = pars["tau"];
    double alpha    = pars["alpha"];
    arma::vec v     = pars["v"];
    arma::uvec z    = pars["z"];
    // Unexchangeable params
    arma::mat beta_unexch = pars["beta_unexch"];
    arma::vec tau_unexch  = pars["tau_unexch"];
    arma::uvec eps0 = pars["eps0"];
    arma::uvec z0 = pars["z0"];
    arma::vec v0 = pars["v_unexch"];
    double alpha_unexch = pars["alpha_unexch"];
    double pexch = pars["pexch"];
    Pars params(
        data, histdata
    , beta, tau, alpha, v, z
    , beta_unexch, tau_unexch, alpha_unexch, v0, pexch, eps0, z0
    );
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
    // Base measure for unexchangeable parameters
    arma::vec beta_unexch_mean = basemeasure["beta_unexch_mean"];
    arma::mat beta_unexch_prec = basemeasure["beta_unexch_prec"];
    double tau_unexch_shape = basemeasure["tau_unexch_shape"];
    double tau_unexch_rate  = basemeasure["tau_unexch_rate"];
    double alpha_unexch_shape = basemeasure["alpha_unexch_shape"];
    double alpha_unexch_rate = basemeasure["alpha_unexch_rate"];
    double alpha_unexch_lower = basemeasure["alpha_unexch_lower"];
    double alpha_unexch_upper = basemeasure["alpha_unexch_upper"];
    double pexch_shape1     = basemeasure["pexch_shape1"];
    double pexch_shape2     = basemeasure["pexch_shape2"];
    double pexch_lower      = basemeasure["pexch_lower"];
    double pexch_upper      = basemeasure["pexch_upper"];
    double pexch_priorprob  = basemeasure["pexch_priorprob"];
    BaseMeasure bm(
        beta_mean, beta_prec, tau_shape, tau_rate
    , alpha_shape, alpha_rate, alpha_lower, alpha_upper 
    , beta_unexch_mean, beta_unexch_prec, tau_unexch_shape, tau_unexch_rate
    , alpha_unexch_shape, alpha_unexch_rate, alpha_unexch_lower, alpha_unexch_upper 
    , pexch_shape1, pexch_shape2, pexch_lower, pexch_upper, pexch_priorprob
    );
    return bm;
  }
} // namespace TTE_BNPDB

#endif // TTE_BNPDB_H
