#include "helperfuns.h"
#include <progress.hpp>
#include "tteData.h"
#include "tte_bnpdb.h"


//' C++ sampler of BNPDB method for time-to-event data
//' 
//' @param data_list a `list` with the elements (y, event, X)c orresponding to the current data
//' @param histdata_list a `list` with the elements (y, event, X) corresponding to the external data
//' @param basemeasure_list a `list` giving the base measure and priors
//' @param inits_list a `list` giving initial values
//' @param niter number of MCMC samples post burn-in and thinning
//' @param nburnin burn-in size
//' @param nthin thin size
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
  tteHistData histdata         = list2tteHistData(histdata_list);
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
