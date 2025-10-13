#include "helperfuns.h"

//' Numerically stable log(sum(exp(x)))
//' @param x arma::vector
//' @return scalar giving log(sum(exp(x)))
//' @name logsumexp
//' @noRd
double logsumexp(arma::vec const& x) {
  double M = x.max();
  return M + log(sum(exp(x - M)));
}


//'Categorical sampler
//' @param log probabilities for categorical distribution
//' @return integer giving a random draw from a categorical distribution
//' @name rcat
//' @noRd
 int rcat(const arma::vec& probs) {
   // Generate a cumulative sum of the probabilities
   arma::vec cumulative_probs = arma::cumsum(probs);
   
   // Generate a random number between 0 and 1
   double u = arma::randu();
   
   // Find the first index where the cumulative probability exceeds the random number
   return arma::as_scalar(arma::find(cumulative_probs > u, 1)); 
 } 

 


//' stick breaking weights from stick breaking variables
//' @param x:   K-1 dim arma::vec of stick-breaking [Beta(1, alpha)] variables
//' @param log: if true, returns result on log scale
//' 
//' @return K dim arma::vec of log of stick breaking weights
//' @noRd
 arma::vec stick_breaking(const arma::vec& x, bool log = false) {
   int K            = x.n_elem + 1;
   arma::vec logx(K, arma::fill::zeros);
   arma::vec log1mx(K, arma::fill::zeros);
   
   logx.head(K - 1)   = arma::log(x);
   log1mx.tail(K - 1) = arma::log1p(-x);
   
   arma::vec logw = logx + arma::cumsum(log1mx);
   
   return log ? logw : arma::exp(logw);
 }   


 
//' Truncated gamma random variable generation
//' 
//' Generates truncated gamma random variates via the inverse-CDF method
//' 
//' @param shape shape parameter for gamma distribution
//' @param rate rate parameter for gamma distribution
//' @param a lower truncation point (non-negative integer)
//' @param b upper truncation point (non-negative integer; b > a)
//' @return single draw from truncated gamma distribution
//' @noRd
 double rtgamma_cpp(double shape, double rate, double a, double b) {
   double scale = 1. / rate;
   
   // If untruncated, use rgamma sampler
   if ( a <= 0 && b == R_PosInf )
     return R::rgamma(shape, scale);
   
   // If truncated only from below and shape is > 1, use sampler used in BART:
   // Gentle J. (2013) Random number generation and Monte Carlo methods. Springer, New York, NY.
   else if ( shape > 1. && a > 0 && b == R_PosInf ) {
     
     // Based off of BART R package
     double a_scale = a * rate;
     double shape_shift = shape - 1.;
     double lambda = 0.5 * (a_scale - shape + sqrt(pow(a_scale - shape, 2.) + 4. * a_scale)) / a_scale;
     double lambda_shift = 1. - lambda;
     double C = 1. + log(lambda_shift / shape_shift);
     
     double y = 0.0;
     double x = 0.0;
     double c = 1.0;
     
     // Sampling loop
     while (c > x) { // Do at least once
       x = R::rexp(1.);
       y = a_scale + R::rexp(1.) / lambda;
       c = lambda_shift * y - shape_shift * (log(y) + C);
     }
     return y / rate;
   } 
   
   // Otherwise, perform slice sampling
   else {
     double x = a + 0.1;
     double w = 0.5;   // Step size for slice sampling, can be adjusted
     double f_x, u, L, R, v, x_new;
     
     // Step 1: Draw a vertical level uniformly from (0, f(x))
     f_x = R::dgamma(x, shape, 1.0 / rate, false);  // Gamma density at x
     u = R::runif(0, f_x);
     
     // Step 2: Create a horizontal interval (L, R) around x
     L = x - w * R::runif(0, 1);
     R = x + w * R::runif(0, 1);
     
     // Truncate the interval to be within (a, b)
     L = std::max(a, L);
     R = std::min(b, R);
     
     // Step 3: Sample from the interval while f(x) < u
     do {
       v = R::runif(0, 1);
       x_new = L + v * (R - L);
       
       // Adjust the interval
       if (x_new < x) {
         L = x_new;
       } else {
         R = x_new;
       }
     } while (R::dgamma(x_new, shape, scale, false) < u);
     return x_new;
   }
   return -1; // never reached
 }  

//' Generate from Dirichlet(1, ..., 1)
//' @param n number of parameters
//' @param log whehter to return result on log scale
//' 
//' @return Dirichlet(1, ..., 1) vector on specified scale
//' @noRd 
 arma::vec rdirichlet_rng(int const& n, bool const& log = false) {
   arma::vec log_gamma(n, arma::fill::zeros);
   for ( int i = 0; i < n; i++ )
     log_gamma(i) = std::log(R::rgamma(1, 1));
   log_gamma -= logsumexp(log_gamma);
   return log ? log_gamma : exp(log_gamma);
 } 
 
 
 
 
 //' Compute conditional survival function based on DPMM regression
 //'
 //' @param Mu nxK matrix giving mean of each observation for each class
 //' @param sigma K-dim vector giving standard deviation for each class
 //' @param logw K-dim vector giving log weight for each class
 //' @param logtimes T-dim vector specifying where to report the survival function
 //'
 //' @return T x n matrix of log survival probabilities
 //' @noRd
 arma::mat compute_survival_conditional(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logtimes, bool log
 ) {
   int T = logtimes.n_elem;
   int n = Mu.n_rows;
   int K = Mu.n_cols;
   arma::vec log_survival_given_xi(K);
   arma::vec log_survival_given_x(n);
   arma::mat out(T, n);
   
   for (int j = 0; j < T; j++) {
     for (int i = 0; i < n; i++) {
       for (int k = 0; k < K; k++) {
         log_survival_given_xi[k] = logw[k] + R::pnorm(logtimes[j], Mu(i,k), sigma[k], 0, 1);
       }
       log_survival_given_x[i] = logsumexp(log_survival_given_xi);
     }
     out.row(j) = log_survival_given_x.t();
   }
   if (log)
     return out;
   return arma::exp(out);
 }

//' Compute marginal survival function based on DPMM regression
//'
//' @param Mu nxK matrix giving mean of each observation for each class
//' @param sigma K-dim vector giving standard deviation for each class
//' @param logw K-dim vector giving log weight for each class
//' @param logomega n-dim vector giving log of Dirichlet random variable
//' @param logtimes T-dim vector specifying where to report the survival function
//' @param log if true, reports log survival function
//'
//' @return T-dim vector giving marginal survival function on specified scale
//' @noRd
 arma::vec compute_survival_bb(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logomega, arma::vec const& logtimes, bool log
 ) {
   arma::mat log_survival_x = compute_survival_conditional(Mu, sigma, logw, logtimes, true);
   int T = logtimes.n_elem;
   arma::vec log_survival(T);
   for (int j = 0; j < T; j++) {
     log_survival[j] = logsumexp(logomega + log_survival_x.row(j).t());
   }
   if (log) return log_survival;
   return arma::exp(log_survival);
 } 

//' Compute conditional RMST function based on DPMM regression
//'
//' @param Mu nxK matrix giving mean of each observation for each class
//' @param sigma K-dim vector giving standard deviation for each class
//' @param logw K-dim vector giving log weight for each class
//' @param logtimes T-dim vector specifying which logtimes to evaluate the RMST
//'
//' @return T x n matrix of log RMST values
//' @noRd
 arma::mat compute_rmst_conditional(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logtimes, bool log
 ) {
   int T = logtimes.n_elem;
   int n = Mu.n_rows;
   int K = Mu.n_cols;
   arma::vec log_rmst_given_xi(K);
   arma::vec log_rmst_given_x(n);
   arma::mat log_rmst_x(T, n);
   arma::vec log_rmst_comp(2);
   double logt_std;
   
   for (int j = 0; j < T; j++) {
     for (int i = 0; i < n; i++) {
       for (int k = 0; k < K; k++) {
         logt_std = (logtimes[j] - Mu(i,k)) / sigma[k];
         log_rmst_comp[0] = logtimes[j] + R::pnorm(logt_std, 0.0, 1.0, 0, 1);
         log_rmst_comp[1] = Mu(i,k) + 0.5 * sigma[k] * sigma[k] + R::pnorm(logt_std - sigma[k], 0.0, 1.0, 1, 1);
         log_rmst_given_xi[k] = logw[k] + logsumexp(log_rmst_comp);
       }
       log_rmst_given_x[i] = logsumexp(log_rmst_given_xi);
     }
     log_rmst_x.row(j) = log_rmst_given_x.t();
   }
   if ( log ) return log_rmst_x;
   else return arma::exp(log_rmst_x);
 } 

//' Compute marginal RMST function based on DPMM regression
//'
//' @param Mu nxK matrix giving mean of each observation for each class
//' @param sigma K-dim vector giving standard deviation for each class
//' @param logw K-dim vector giving log weight for each class
//' @param logomega n-dim vector giving log of Dirichlet random variable
//' @param logtimes T-dim vector specifying which logtimes to evaluate the RMST
//' @param log if true, reports log RMST
//'
//' @return T-dim vector giving marginal RMST on specified scale
//' @noRd
 arma::vec compute_rmst_bb(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logomega, arma::vec const& logtimes, bool log
 ) {
   arma::mat log_rmst_x = compute_rmst_conditional(Mu, sigma, logw, logtimes, true);
   int T = logtimes.n_elem;
   arma::vec log_rmst(T);
   for (int j = 0; j < T; j++) {
     log_rmst[j] = logsumexp(logomega + log_rmst_x.row(j).t());
   }
   if (log) return log_rmst;
   return arma::exp(log_rmst);
 } 

//' Compute conditional hazard function based on DPMM regression
//' 
//' @param Mu nxK matrix giving mean of each observation for each class
//' @param sigma K-dim vector giving standard deviation for each class
//' @param logw K-dim vector giving log weight for each class
//' @param logtimes T-dim vector specifying which logtimes to evaluate the hazard
//'
//' @return T x n matrix of log hazard values
//' @noRd
 arma::mat compute_hazard_conditional(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logtimes, bool log
 ) {
   int T = logtimes.n_elem;
   int n = Mu.n_rows;
   int K = Mu.n_cols;
   
   arma::vec log_density_given_xi(K);
   arma::vec log_survival_given_xi(K);
   arma::vec log_density_given_x(n);
   arma::vec log_survival_given_x(n);
   arma::mat log_hazard_x(T, n);
   
   for (int j = 0; j < T; j++) {
     for (int i = 0; i < n; i++) {
       for (int k = 0; k < K; k++) {
         log_density_given_xi[k]  = logw[k] - logtimes[j] + R::dnorm(logtimes[j], Mu(i,k), sigma[k], 1);
         log_survival_given_xi[k] = logw[k] + R::pnorm(logtimes[j], Mu(i,k), sigma[k], 0, 1);
       }
       log_density_given_x[i]  = logsumexp(log_density_given_xi);
       log_survival_given_x[i] = logsumexp(log_survival_given_xi);
     }
     log_hazard_x.row(j) = (log_density_given_x - log_survival_given_x).t();
   }
   if ( log ) return log_hazard_x;
   else return arma::exp(log_hazard_x);
 }

//' Compute marginal hazard function based on DPMM regression
//'
//' @param Mu nxK matrix giving mean of each observation for each class
//' @param sigma K-dim vector giving standard deviation for each class
//' @param logw K-dim vector giving log weight for each class
//' @param logomega n-dim vector giving log of Dirichlet random variable
//' @param logtimes T-dim vector specifying which logtimes to evaluate the hazard
//' @param log if true, reports log hazard
//'
//' @return T-dim vector giving marginal hazard function on specified scale
//' @noRd
 arma::vec compute_hazard_bb(
     arma::mat const& Mu, arma::rowvec const& sigma, arma::rowvec const& logw,
     arma::vec const& logomega, arma::vec const& logtimes, bool log
 ) {
   arma::mat log_hazard_x = compute_hazard_conditional(Mu, sigma, logw, logtimes, log = true);
   int T = logtimes.n_elem;
   arma::vec log_hazard(T);
   for (int j = 0; j < T; j++) {
     log_hazard[j] = logsumexp(logomega + log_hazard_x.row(j).t());
   }
   if (log) return log_hazard;
   return arma::exp(log_hazard);
 } 

 
 //' Convert Rcpp list of data to a tteData object
 //' @keywords internal
 //' @noRd
 tteData list2tteData(Rcpp::List const& data) {
   arma::vec y      = data["y"];
   arma::uvec event = data["event"];
   arma::mat X      = data["X"];
   tteData dat(y, event, X);
   return dat;
 }


//' Convert Rcpp list of historical data to a tteHistData object
//' @keywords internal
//' @noRd
 tteHistData list2tteHistData(Rcpp::List const& histdata) {
   arma::vec y       = histdata["y"];
   arma::uvec event  = histdata["event"];
   arma::mat X       = histdata["X"];
   arma::mat Xunexch = histdata["Xunexch"];
   tteHistData histdat(y, event, X, Xunexch);
   return histdat;
 } 



//' Sample from truncated beta distribution
//' @keywords internal
//' @noRd
//' @export
 double rtbeta(double shape1, double shape2, double lower, double upper) {
   lower = std::max(0.0, lower);
   upper = std::min(1.0, upper);
   
   if (lower == upper)
     return lower;
   else if (lower == 0 && upper == 1)
     return R::rbeta(shape1, shape2);
   else if (shape1 == 1 && shape2 == 1)
     return R::runif( lower, upper );
   else {
     // Use inverse transform
     double Fupper = R::pbeta(upper, shape1, shape2, 1, 0);
     double Flower = R::pbeta(lower, shape1, shape2, 1, 0);
     double u      = R::runif(Flower, Fupper);
     return R::qbeta(u, shape1, shape2, 1, 0);
   }
   return(-1.0);   // never reached
 }



//' Numerically stable log(exp(logb) - exp(loga))
//' @param logb larger logged number
//' @param loga smaller logged number
//' @noRd
double logdiffexp(double logb, double loga) {
 if (loga == R_NegInf) return logb;
 if (loga > logb) return R_NaN;       // invalid: negative difference
 if (loga == logb) return R_NegInf;   // log(0)
 
 double d = loga - logb;
 // exp(d) underflows safely for very negative d
 return logb + log1p(-std::exp(d));   // numerically stable
}

// [[Rcpp::export]]
double sample_tbeta_mixture_posterior(
    double shape1, double shape2, double lower, double upper,
    double prior_prob, int nexch, int nunexch
) {
  // Return point mass if prior_prob == 1
  if (prior_prob == 1) return upper;
  
  // Calculate posterior parameters
  double shape1_post = shape1 + nexch;
  double shape2_post = shape2 + nunexch;
  
  // If prior_prob is zero or upper bound is 1 and there are unexchangeable
  // subjects, sample from truncated beta ( it is impossible for pexch to be 
  // equal to 1  if there are unexchangeable subjects)
  if (prior_prob == 0 || (nunexch > 0 && upper == 1) )
    return rtbeta(shape1_post, shape2_post, lower, upper);
  
  // Otherwise, we have a mixture of a (possibly truncated) beta distribution
  // and a point mass at the upper bound (<=1)
  
  // The following are helpful to compute the updated mixture weight
  double lognc_prior = R::lbeta(shape1, shape2);
  double lognc_post  = R::lbeta(shape1_post, shape2_post);
  double logD_prior, logD_post;
  
  if (lower == 0 && upper == 1) {
    logD_prior = 0;
    logD_post  = 0;
  } else if (lower == 0 && upper < 1) {
    logD_prior = R::pbeta(upper, shape1, shape2, 1, 1);
    logD_post  = R::pbeta(upper, shape1_post, shape2_post, 1, 1);
  } else if (lower > 0 && upper == 1) { 
    logD_prior = R::pbeta(lower, shape1, shape2, 0, 1);  // log[1 - F(a)]
    logD_post  = R::pbeta(lower, shape1_post, shape2_post, 0, 1);
  } else { 
    logD_prior = logdiffexp(R::pbeta(upper, shape1, shape2, 1, 1),
                            R::pbeta(lower, shape1, shape2, 1, 1));
    logD_post  = logdiffexp(R::pbeta(upper, shape1_post, shape2_post, 1, 1),
                            R::pbeta(lower, shape1_post, shape2_post, 1, 1));
  }
  lognc_prior += logD_prior;
  lognc_post  += logD_post;
  
  // Compute log components (unnormalized) of each mixture weight. The first
  // component refers to the point mass at `upper` and the second component
  // refers to the truncated beta
  arma::vec logpp_comp(2);
  logpp_comp[0] = log(prior_prob);
  if (upper < 1)
    logpp_comp[0] += nexch * log(upper) + nunexch * log1p(-upper);
  logpp_comp[1] = log1p(-prior_prob) + lognc_post - lognc_prior;
  
  // Compute mixture probability; draw a random binary variable z; then draw
  // from point mass ( if z == 1 ) or truncated beta posterior ( if z == 0 )
  double p = exp(logpp_comp[0] - logsumexp(logpp_comp));
  int z = R::rbinom(1, p);
  return (z == 1) ? upper : rtbeta(shape1_post, shape2_post, lower, upper);
} 




