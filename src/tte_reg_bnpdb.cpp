#include "helperfuns.h"
#include <progress.hpp>
#include "tteData.h"
#include "tte_bnpdb.h"


//' C++ sampler of DPMM method for time-to-event data
//' 
//' @param y log of observed event time for current data
//' @param X design matrix for current data
//' @param event event indicator for current data (1 = event; 0 = right censored)
//' @param K number of stick breaking components for exchangeable density
//' @param beta_mean mean hyperparameter for normal base measure on regression coefficients conditional on precision for exchangeable density
//' @param beta_cov scale hyperparameter for normal base measure on regression coefficients conditional on precision for exchangeable density
//' @param tau_shape shape hyperparameter for gamma base measure on precision parameter for exchangeable density
//' @param tau_rate rate hyperparameter for gamma base measure on precision parameter for exchangeable density
//' @param alpha_shape shape hyperparameter for truncated gamma prior on DP concentration parameter for exchangeable density
//' @param alpha_rate rate hyperparameter for truncated gamma prior on DP concentration parameter for exchangeable density
//' @param alpha_lower lower bound for truncated gamma prior on DP concentration parameter for exchangeable density
//' @param alpha_upper upper bound for truncated gamma prior on DP concentration parameter for exchangeable density
//' @param niter number of desired posterior samples post burn-in
//' @param nburnin number of burn-in samples
//' @param thin number of samples to thin post burn-in
//' 
//' @return a list giving posterior samples
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List tte_reg_bnpdb_cpp(
    Rcpp::List const& data_list, Rcpp::List const& histdata_list
    , Rcpp::List const& basemeasure_list, Rcpp::List const& inits_list
    , int const& niter, int const& nburnin, int const& nthin
) {
  // Initialize data, historical data base measure, and parameters (pars) objects
  tteData data                 = list2tteData(data_list);
  tteData histdata             = list2tteData(histdata_list);
  TTE_BNPDB::BaseMeasure bm    = TTE_BNPDB::list2basemeasure(basemeasure_list);
  TTE_BNPDB::Pars pars         = TTE_BNPDB::list2pars(data, histdata, inits_list);
  
  // Initialize container for results
  arma::cube beta(data.p, pars.K, niter, arma::fill::zeros);
  arma::mat  sigma(pars.K, niter, arma::fill::zeros);
  arma::vec  alpha(niter, arma::fill::zeros);
  arma::mat  w(pars.K, niter, arma::fill::zeros);
  arma::umat eps0(histdata.n, niter, arma::fill::zeros);
  arma::vec  pexch(niter, arma::fill::zeros);
  
  // Initialize progress bar
  Progress progress(nburnin + niter, true);
  
  // Conduct MCMC
  for ( int i = 0; i < nburnin; i++ ) {
    // Check for user interrupt
    if (Progress::check_abort()) return Rcpp::List::create();
    
    // Update MCMC parameters
    pars.update_pexch(bm);
    pars.update_concentration_parameter_exch(bm);
    pars.update_concentration_parameter_unexch(bm);
    pars.update_sb_variables_exch();
    pars.update_sb_variables_unexch();
    pars.update_classes_current(data);
    pars.update_classes_historical(histdata);
    pars.update_regression_params_exch(data, histdata, bm);
    pars.update_regression_params_unexch(histdata, bm);
    pars.update_censored_vars_current(data);
    pars.update_censored_vars_historical(histdata);
    progress.increment();
  } 
  
  for ( int i = 0; i < niter; i++ ) {
    // Check for user interrupt
    if (Progress::check_abort()) return Rcpp::List::create();
    
    // Keep only the thinned samples
    for ( int j = 0; j < nthin; j++ ) {
      // Update MCMC parameters
      pars.update_pexch(bm);
      pars.update_concentration_parameter_exch(bm);
      pars.update_concentration_parameter_unexch(bm);
      pars.update_sb_variables_exch();
      pars.update_sb_variables_unexch();
      pars.update_classes_current(data);
      pars.update_classes_historical(histdata);
      pars.update_regression_params_exch(data, histdata, bm);
      pars.update_regression_params_unexch(histdata, bm);
      pars.update_censored_vars_current(data);
      pars.update_censored_vars_historical(histdata);
    } 
    // Store results
    beta.slice(i) = pars.beta;
    sigma.col(i)  = pars.sigma;
    alpha(i)      = pars.alpha;
    w.col(i)      = exp(pars.logw);
    eps0.col(i)   = pars.eps0;
    pexch(i)      = pars.pexch;
    
    progress.increment(); // Update progress bar
  } 
  
  sigma    = sigma.t();
  w        = w.t();
  
  Rcpp::List res = Rcpp::List::create(
    Named("beta")     = beta
  , Named("sigma")    = sigma
  , Named("alpha")    = alpha
  , Named("w")        = w
  , Named("pexch")    = pexch
  , Named("eps0")     = eps0.t()
  );
  return res;
} 
