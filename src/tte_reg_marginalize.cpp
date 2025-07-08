#include "helperfuns.h"
#include <progress.hpp>

//' Obtain marginal quantities when design matrix consists only of strata.
//' 
//' 
//' @param X design matrix with S rows and S columns, S = number of strata.
//' @param samples a named `list` giving posterior samples of the DPMM
//' @param log_survival_times log of times to evaluate survival function
//' @param log_survival boolean; if true: outputs log of survival function
//' @param log_hazard_times log of time points to evaluate hazard function
//' @param log_hazard boolean; if true: outputs log of hazard function
//' @param log_rmst_times log of times to evaluate RMST function
//' @param log_rmst boolean; if true: outputs log of RMST function
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
 Rcpp::List marginalize_strataonly(
     arma::mat const& X, Rcpp::List const& samples
    , arma::vec const& log_survival_times, bool log_survival
    , arma::vec const& log_hazard_times, bool log_hazard
    , arma::vec const& log_rmst_times, bool log_rmst
 ) {
   arma::cube beta  = Rcpp::as<arma::cube>(samples["beta"]);
   arma::mat  sigma = Rcpp::as<arma::mat>(samples["sigma"]);
   arma::mat  logw  = Rcpp::as<arma::mat>(samples["w"]);
   logw             = arma::log(logw);
   
   int M      = logw.n_rows;                // number of MCMC samples
   int m_surv = log_survival_times.n_elem;  // number of survival times
   int m_rmst = log_rmst_times.n_elem;      // number of RMST times
   int m_haz  = log_hazard_times.n_elem;    // number of hazard times
   int S      = X.n_rows;                   // number of strata (>= 1)
   
   // Return empty list if nothing to evaluate
   if (m_surv == 0 && m_rmst == 0 && m_haz == 0)
     return Rcpp::List::create();
   
   // Initialize output matrices (time points × MCMC samples)
   arma::cube survival_arr(M, m_surv, S, arma::fill::zeros);
   arma::cube rmst_arr(M, m_rmst, S, arma::fill::zeros);
   arma::cube hazard_arr(M, m_haz, S, arma::fill::zeros);
   
   // Initialize progress bar
   Progress progress(M, true);
   
   // Loop thorugh MCMC samples
   for ( int m = 0; m < M; m++ ) {
     
     // Check if user interruption
     if (Progress::check_abort()) return Rcpp::List::create();
     
     // Compute stratum-specific means
     arma::mat Mu = X * beta.slice(m);  // S x K matrix
     
     if ( m_surv > 0 ) {
       arma::mat surv = compute_survival_conditional(Mu, sigma.row(m), logw.row(m), log_survival_times, log_survival);
       survival_arr.row(m) = surv;
     }
     
     if ( m_rmst > 0 ) {
       arma::mat rmst = compute_rmst_conditional(Mu, sigma.row(m), logw.row(m), log_rmst_times, log_rmst);
       rmst_arr.row(m) = rmst;
     }
     
     if ( m_haz > 0 ) {
       arma::mat haz     = compute_hazard_conditional(Mu, sigma.row(m), logw.row(m), log_hazard_times, log_hazard);  // T x S matrix
       hazard_arr.row(m) = haz;
     }
     
     progress.increment();
     
   }
   
   return Rcpp::List::create(
     Rcpp::Named("survival") = survival_arr,
     Rcpp::Named("rmst")     = rmst_arr,
     Rcpp::Named("hazard")   = hazard_arr
   );
 }


//' Marginalize over strata for time-to-event regression using posterior samples
//'
//' @param Xlist a list of S design matrices (one per stratum), each with n rows
//' @param samples_list a list of length S, each element a named list of posterior samples (`beta`, `sigma`, `w`)
//' @param log_survival_times log of times to evaluate survival function
//' @param log_survival boolean; if true: outputs log of survival function
//' @param log_hazard_times log of time points to evaluate hazard function
//' @param log_hazard boolean; if true: outputs log of hazard function
//' @param log_rmst_times log of times to evaluate RMST function
//' @param log_rmst boolean; if true: outputs log of RMST function
//'
//' @return a named list of three 3D arrays (`survival`, `rmst`, `hazard`) with dimensions (time × MCMC × S)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
 Rcpp::List marginalize_tte_gcomp(
     Rcpp::List const& Xlist,
     Rcpp::List const& samples_list,
     arma::vec const& log_survival_times,
     bool log_survival,
     arma::vec const& log_hazard_times,
     bool log_hazard,
     arma::vec const& log_rmst_times,
     bool log_rmst
 ) {
   int S = Xlist.size();  // Number of strata
   if (S == 0) return Rcpp::List::create();
   
   // Extract MCMC dimensions from first stratum
   Rcpp::List samples0 = samples_list[0];
   arma::cube beta0    = Rcpp::as<arma::cube>(samples0["beta"]);
   int M = beta0.n_slices;
   int n =  Rcpp::as<arma::mat>(Xlist[0]).n_rows;
   
   int m_surv = log_survival_times.n_elem;
   int m_rmst = log_rmst_times.n_elem;
   int m_haz  = log_hazard_times.n_elem;
   
   if (m_surv == 0 && m_rmst == 0 && m_haz == 0)
     return Rcpp::List::create();
   
   arma::cube survival_arr(M, m_surv, S, arma::fill::zeros);
   arma::cube rmst_arr(M, m_rmst, S, arma::fill::zeros);
   arma::cube hazard_arr(M, m_haz, S, arma::fill::zeros);
   
   Progress progress(M*S, true);
   
   // Pre-generate the log Dirichlet random variables to allow looping
   // over strata first
   arma::mat logomega(n, M, arma::fill::zeros);
   for ( int m = 0; m < M; m++ ) {
     logomega.col(m) = rdirichlet_rng(n, true); 
   }
   
   // Loop over strata
   for (int s = 0; s < S; s++) {
     // Extract quantities for stratum s
     Rcpp::List samples = samples_list[s];
     arma::cube beta    = Rcpp::as<arma::cube>(samples["beta"]);
     arma::mat sigma    = Rcpp::as<arma::mat>(samples["sigma"]);
     arma::mat logw     = arma::log(Rcpp::as<arma::mat>(samples["w"]));
     arma::mat X        = Rcpp::as<arma::mat>(Xlist[s]);
     
     // loop over MCMC samples
     for (int m = 0; m < M; m++) {
       if (Progress::check_abort()) return Rcpp::List::create();
       
       arma::mat Mu = X * beta.slice(m);  // n × K matrix of stratum means
       
       if (m_surv > 0) {
         arma::vec surv = compute_survival_bb(Mu, sigma.row(m), logw.row(m), logomega.col(m), log_survival_times, log_survival);
         survival_arr.slice(s).row(m) = surv.t();
       } 
       
       if (m_rmst > 0) {
         arma::vec rmst = compute_rmst_bb(Mu, sigma.row(m), logw.row(m), logomega.col(m), log_rmst_times, log_rmst);
         rmst_arr.slice(s).row(m) = rmst.t();
       } 
       
       if (m_haz > 0) {
         arma::vec haz = compute_hazard_bb(Mu, sigma.row(m), logw.row(m), logomega.col(m), log_hazard_times, log_hazard);
         hazard_arr.slice(s).row(m) = haz.t();
       }
       
       progress.increment();
       
     } 
   }
   
   return Rcpp::List::create(
     Rcpp::Named("survival") = survival_arr,
     Rcpp::Named("rmst")     = rmst_arr,
     Rcpp::Named("hazard")   = hazard_arr
   );
 } 


