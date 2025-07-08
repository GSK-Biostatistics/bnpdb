
##--------------------------------------------------------------------------
## This file contains functions to check the priors / base measures for 
## continuous DPMM models, using data-driven ones if the user left them
## unspecified. 
##
## This file contains helper functions, and is not part of the main output
## of bnpdb.
##--------------------------------------------------------------------------


#' Find gamma distribution that minimizes KL divergence based on log-normal
#' distribution
#' 
#' @import stats
#' @keywords internal
#' @noRd
lognormal2gamma <- function(meanlog, sdlog) {
  # Log-normal density
  dln <- function(x, log = FALSE) {
    dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = log)
  }
  
  # KL divergence from LogNormal to Gamma
  kl_div <- function(par) {
    shape <- par[1]
    rate  <- par[2]
    if (shape <= 0 || rate <= 0) return(Inf)
    
    integrand <- function(x) {
      res         <- numeric(length(x))
      indx0       <- x <= 0
      res[indx0]  <- 0
      x           <- x[!indx0]
      
      log_dln_x   <- dln(x, TRUE)
      log_ratio   <- log_dln_x - dgamma(x, shape, rate, log = TRUE)
      res[!indx0] <- exp(log_dln_x + log_ratio)
      
      return(res)
    }
    
    lower <- qlnorm(0.01, meanlog, sdlog)
    upper <- qlnorm(0.99, meanlog, sdlog)
    result <- tryCatch({
      integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)$value
    }, error = function(e) 1e6)
    
    return(result)
  }
  
  # Initial guess using method of moments
  mean_ln    <- exp(meanlog + 0.5 * sdlog^2)
  var_ln     <- (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
  shape_init <- mean_ln^2 / var_ln
  rate_init  <- mean_ln / var_ln
  par_init   <- c(shape_init, rate_init)
  
  # Optimization
  opt <- optim(par_init, kl_div, method = "L-BFGS-B",
               lower = c(1e-4, 1e-4))
  if ( opt$convergence != 0 )
    warning("Could not minimize KL divergence")
  
  res = list(shape = opt$par[1], rate = opt$par[2], kl = opt$value)
  attr(res, 'optim') = opt
  res
}

check_prior_list_continuous <- function(
    fit, alpha_prior, beta_bm, tau_bm, X, external = FALSE
) {
  p = ncol(X)
  n = nrow(X)
  hyperparams = list()
  ## Hyperparameters for alpha
  if ( family(alpha_prior) == 'truncated' ) {
    hyperparams$alpha_lower = parameters( alpha_prior )[['lower']]
    hyperparams$alpha_upper = parameters( alpha_prior )[['upper']]
    dist = parameters(alpha_prior)[['dist']]
    stopifnot("If alpha_prior is truncated, it must be a truncated gamma distribution" = (family(dist) == 'gamma') )
    hyperparams$alpha_shape = parameters(dist)[['shape']]
    hyperparams$alpha_rate  = parameters(dist)[['rate']]
  } else if ( family(alpha_prior) == 'gamma' ) {
    hyperparams$alpha_lower = 0
    hyperparams$alpha_upper = Inf
    hyperparams$alpha_shape = parameters(alpha_prior)[['shape']]
    hyperparams$alpha_rate  = parameters(alpha_prior)[['rate']]
  } else stop("alpha_prior must be a gamma prior or a truncated gamma prior.")
  
  ##------------------------------
  ## Hyperparameters for beta
  ##------------------------------
  
  ## If base measure for beta is not specified, elicit a multivariate normal base
  ## measure with mean = MLE and covariance = 5 events worth of information
  if ( is.null(beta_bm) ) {
    nevents   = sum(fit$y[, 2])
    beta_mean = as.numeric(coef(fit))
    names(beta_mean) = NULL
    beta_cov  = nevents / 5 * vcov(fit)[1:p, 1:p]
    beta_bm   = dist_multivariate_normal(mu = list(beta_mean), sigma = list(beta_cov))
  }
  stopifnot( 
    "beta_bm not of correct length. Try specifying beta_bm = dist_normal(mu = list(<vector of means>), sigma = list(<covariance matrix>)) where each list contains one element."
    = (length(beta_bm) == 1)
  )
  stopifnot( "Prior on beta must be multivariate normal" = (family(beta_bm) == 'mvnorm') )
  
  hyperparams$beta_mean = mean(beta_bm)[1, ]
  hyperparams$beta_cov  = covariance(beta_bm)[[1]]
  hyperparams$beta_prec = solve(hyperparams$beta_cov)
  stopifnot( "Length of prior mean must match number of columns of design matrix" = length(hyperparams$beta_mean) == p )
  stopifnot( "Prior covariance must be square with dimensions equal to the number of columns of the design matrix" = all( nrow(hyperparams$beta_cov) == p, ncol(hyperparams$beta_cov) == p ) )
  stopifnot( "Prior covariance matrix for beta must be positive definite" = all(eigen(hyperparams$beta_cov)$values > 0) )
  
  ##------------------------------
  ## Hyperparameters for tau
  ##------------------------------
  
  ## If base measure for tau is not specified, elicit a gamma distribution that
  ## minimizes the KL divergence between the log-normal distribution implied
  ## by the MLE for the standard deviation; using the fact that
  ##  sigma ~ LN(mu, sigma) => tau := 1 / sigma^2 ~ LN(-2 * mu, 2 * sigma)
  
  if ( is.null ( tau_bm ) ) {
    nevents      = sum(fit$y[, 2])
    tau_ln_mean  = -2 * log( fit$scale )
    tau_ln_sd    = 2 * sqrt(nevents / 5 * vcov(fit)['Log(scale)', 'Log(scale)'])
    ln2gamma     = lognormal2gamma(tau_ln_mean, tau_ln_sd)
    tau_bm       = dist_gamma(ln2gamma$shape, ln2gamma$rate)
  }
  
  stopifnot("Base measure for tau must be a gamma distribution" = family(tau_bm) == 'gamma')
  
  hyperparams$tau_shape = parameters(tau_bm)[['shape']]
  hyperparams$tau_rate  = parameters(tau_bm)[['rate']]
  
  ## Return list
  res = list(
    hyperparams = hyperparams, alpha_prior = alpha_prior, beta_bm = beta_bm
    , tau_bm = tau_bm
  )
  ## If external flag is TRUE, change the names of hyperparameters and
  ## prior to be of the form x_unexch
  if ( external ) {
    ## Replace hyperparameter names to be x_unexch_...
    names(res$hyperparams) = sub("^([^_]+)_", "\\1_unexch_", names(res$hyperparams))
    ## Replace prior names to be x_unexch_...
    indx = which(names(res) != 'hyperparams')
    names(res)[indx] = sub("^([^_]+)_", "\\1_unexch_", names(res)[indx])
  }
  
  return(res)
  
}


