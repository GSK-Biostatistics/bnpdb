#include "helperfuns.h"
#include <progress.hpp>
#include "tteData.h"
#include "tte_dpmm.h"


//' C++ sampler of DPMM method for time-to-event data
//' 
//' @param niter number of desired posterior samples post burn-in
//' @param nburnin number of burn-in samples
//' @param thin number of samples to thin post burn-in
//' 
//' @return a list giving posterior samples
//' @noRd
// [[Rcpp::export]]
Rcpp::List tte_reg_dpmm_cpp(
    Rcpp::List const& data_list
  , Rcpp::List const& basemeasure_list, Rcpp::List const& inits_list
  , int const& niter, int const& nburnin, int const& nthin
) {
  // Initialize data, base measure, and parameters (pars) objects
  tteData data             = list2tteData(data_list);
  TTE_DPMM::BaseMeasure bm = TTE_DPMM::list2basemeasure(basemeasure_list);
  TTE_DPMM::Pars pars      = TTE_DPMM::list2pars(data, inits_list);
  
  // Initialize container for results
  arma::cube beta(data.p, pars.K, niter, arma::fill::zeros);
  arma::mat  sigma(pars.K, niter, arma::fill::zeros);
  arma::vec  alpha(niter, arma::fill::zeros);
  arma::mat  w(pars.K, niter, arma::fill::zeros);
  
  
  // Initialize progress bar
  Progress progress(nburnin + niter, true);
  
  // Conduct MCMC
  for ( int i = 0; i < nburnin; i++ ) {
    // Check for user interrupt
    if (Progress::check_abort()) return Rcpp::List::create();
    
    // Update MCMC parameters
    pars.update_concentration_parameter(bm);
    pars.update_sb_variables();
    pars.update_classes(data);
    pars.update_regression_params(data, bm);
    pars.update_censored_vars(data);
    progress.increment();
  } 
  
  for ( int i = 0; i < niter; i++ ) {
    // Check for user interrupt
    if (Progress::check_abort()) return Rcpp::List::create();
    
    // Keep only the thinned samples
    for ( int j = 0; j < nthin; j++ ) {
      // Update MCMC parameters
      pars.update_concentration_parameter(bm);
      pars.update_sb_variables();
      pars.update_classes(data);
      pars.update_regression_params(data, bm);
      pars.update_censored_vars(data);
    } 
    // Store results
    beta.slice(i) = pars.beta;
    sigma.col(i)  = pars.sigma;
    alpha(i)      = pars.alpha;
    w.col(i)      = exp(pars.logw);
    
    progress.increment(); // Update progress bar
  } 
  
  sigma    = sigma.t();
  w        = w.t();
  
  Rcpp::List res = Rcpp::List::create(
    Named("beta")     = beta
  , Named("sigma")    = sigma
  , Named("alpha")    = alpha
  , Named("w")        = w
  );
  return res;
} 

