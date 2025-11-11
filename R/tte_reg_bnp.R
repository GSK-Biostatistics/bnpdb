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
    
    lower <- qlnorm(0.001, meanlog, sdlog)
    upper <- qlnorm(0.999, meanlog, sdlog)
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
               lower = c(1.1, 1e-4))
  if ( opt$convergence != 0 )
    warning("Could not minimize KL divergence. Using method of moments")
  
  res = list(shape = opt$par[1], rate = opt$par[2], kl = opt$value)
  attr(res, 'optim') = opt
  res
}



#' Construct data and default base measure to feed onto C++ 
#' 
#' @import stats
#' @import distributional
#' @import survival
#' @keywords internal
#' @noRd
get_bm_tte = function (
    bm = list(), n, mle, vcov, mle0 = NULL, vcov0 = NULL, n0 = NULL, prior_ess = 5
  ) {
  p     = length(mle) - 1
  sigma = vcov * (n / prior_ess)
  if ( is.null(bm$beta_bm) ) {
    bm$beta_bm = distributional::dist_multivariate_normal(
      mu = list( mle[1:p] ), sigma = list( as.matrix(sigma[1:p,1:p]) )
    )
  }
  if ( is.null( bm$tau_bm ) ) {
    ## Fit gamma approximation to log-normal distribution
    logtau_mean  = -2 * mle['Log(scale)']
    logtau_sd    =  2 * sqrt(sigma['Log(scale)', 'Log(scale)'])
    parms        = lognormal2gamma(logtau_mean, logtau_sd)
    bm$tau_bm    = distributional::dist_gamma(parms$shape, parms$rate)
  }
  if ( is.null( bm$alpha_prior ) )
    bm$alpha_prior = distributional::dist_gamma(8, 2)
  
  if ( !is.null(mle0) ) {
    p0 = length(mle0) - 1
    sigma0 = vcov0 * (n0 / prior_ess)
    if ( is.null(bm$beta_unexch_bm) ) {
      bm$beta_unexch_bm = distributional::dist_multivariate_normal(
        mu = list( rep(0, p0) ), sigma = list( as.matrix( sigma0[1:p0,1:p0] )  )
      )
    }
    if ( is.null( bm$tau_unexch_bm ) ) {
      ## Fit gamma approximation to log-normal distribution
      logtau_mean      = -2 * mle0['Log(scale)']
      logtau_sd        =  2 * sqrt(sigma0['Log(scale)', 'Log(scale)'])
      parms            = lognormal2gamma(logtau_mean, logtau_sd)
      bm$tau_unexch_bm = distributional::dist_gamma(parms$shape, parms$rate)
    }
    if ( is.null(bm$alpha_unexch_prior) )
      bm$alpha_unexch_prior = distributional::dist_gamma(8, 2)
    
    if ( is.null(bm$pexch_prior) ) {
      dist1 = distributional::dist_beta(1, 1)
      dist2 = distributional::dist_degenerate(1)
      bm$pexch_prior = distributional::dist_mixture(dist1, dist2, weights = c(0.5, 0.5))
    }
  }
  bm
}



#' Construct default initial values if user does not specify
#' 
#' @import stats
#' @import distributional
#' @import survival
#' @keywords internal
#' @noRd
tte_inits = function(
  inits = list(), bm, cppData, K, cppExtData = NULL, K_unexch = NULL
) {
  X     = cppData$X
  y     = cppData$y
  event = cppData$event
  n     = nrow(X)
  
  ## If beta initial values are not specified, take one component to be the mean
  ## and randomly draw the other components from base measure
  if ( is.null(inits$beta) ) {
    inits$beta = cbind( mean(bm$beta_bm)[1, ], t( generate(bm$beta_bm, K-1)[[1]] ) )
  }
  ## If tau initial values are not specified, take one component to be the mean
  ## and randomly draw the other components from base measure
  if ( is.null(inits$tau) ) {
    inits$tau = c(mean(bm$tau_bm), generate(bm$tau_bm, K-1)[[1]])
  }
  if ( is.null(inits$z) ) {
    mode       = rep(0, n)
    logProbMax = rep(-Inf, n)
    for ( k in 1:K ) {
      mu_k    = X %*% inits$beta[, k]
      sigma_k = 1 / sqrt(inits$tau[k])
      logProb = ifelse(event == 1, dnorm(y, mu_k, sigma_k, log = TRUE), pnorm(y, mu_k, sigma_k, lower.tail = FALSE, log.p = TRUE))
      indx    = logProb > logProbMax
      mode[indx] = k
      logProbMax[indx] = logProb[indx]
    }
    inits$z = mode
  }
  if ( is.null(inits$alpha) )
    inits$alpha = mean(bm$alpha_prior)
  if ( is.null(inits$v) )
    inits$v = rbeta(K-1, 1, inits$alpha)
  
  if ( !is.null(cppExtData) ) {
    y0        = cppExtData$y
    event0    = cppExtData$event
    X0        = cppExtData$X
    X0_unexch = cppExtData$Xunexch
    n0        = length(y0)
    
    if ( is.null(inits$alpha_unexch) ) {
      inits$alpha_unexch = mean(bm$alpha_unexch_prior)
    }
    if ( is.null(inits$v_unexch) ) {
      inits$v_unexch = rbeta(K_unexch-1, 1, inits$alpha_unexch)
    }
    if ( is.null(inits$beta_unexch) ) {
      inits$beta_unexch = cbind( mean(bm$beta_unexch_bm)[1, ], t( generate(bm$beta_unexch_bm, K_unexch-1)[[1]] ) )
    }
    if ( is.null(inits$tau_unexch) ) {
      inits$tau_unexch = c(mean(bm$tau_unexch_bm), generate(bm$tau_unexch_bm, K_unexch-1)[[1]])
    }
    
    ## Obtain log likelihood contributions for historical data to help specify
    ## initial values for `z0` and `eps0`
    if ( is.null(inits$eps0) & !is.null(inits$z0) ) {
      stop("If `z0_unexch` is specified then `eps0` must also be specified")
    } else if ( is.null(inits$eps0) | is.null(inits$z0) ) {
      maxLogLik  = rep(-Inf, n0)
      maxlogLik0 = rep(-Inf, n0)
      z0         = rep(0, n0)
      z0_unexch  = rep(0, n0)
      for ( k in 1:K ) {
        mu_k            = X0 %*% inits$beta[, k]
        sigma_k         = 1 / sqrt(inits$tau[k])
        loglik_k        = ifelse( event0 == 1, dnorm(y0, mu_k, sigma_k, log = TRUE), pnorm(y0, mu_k, sigma_k, log.p = TRUE, lower.tail = FALSE) )
        indx            = loglik_k > maxLogLik
        maxLogLik[indx] = loglik_k[indx]
        z0[indx] = k
      }
      for ( k in 1:K_unexch ) {
        mu_k             = X0_unexch %*% inits$beta_unexch[, k]
        sigma_k          = 1 / sqrt(inits$tau_unexch[k])
        loglik_k         = ifelse( event0 == 1, dnorm(y0, mu_k, sigma_k, log = TRUE), pnorm(y0, mu_k, sigma_k, log.p = TRUE, lower.tail = FALSE) )
        indx             = loglik_k > maxlogLik0
        maxlogLik0[indx] = loglik_k[indx]
        z0_unexch[indx]  = k
      }
      if ( is.null(inits$eps0) )
        inits$eps0      = ifelse(maxLogLik > maxlogLik0, 1, 0)
      inits$z0 = ifelse(inits$eps0 == 1, z0, K + z0_unexch)
    }
    
    if ( is.null(inits$pexch) )
      inits$pexch = mean(bm$pexch_prior)
    
  }
  inits
}

## Function to perform checks for time-to-event BNPDB
tte_check_all_and_export_hyperparams = function(
  cppData, K, bm, inits, cppExtData = NULL, K_unexch = NULL
) {
  p   = ncol(cppData$X)
  n   = nrow(cppData$X)
  hyp = list()
  
  ## Check regression coefficients for current data: base measure and initial values
  stopifnot("`beta_bm` is not specified" = !is.null(bm$beta_bm))
  stopifnot(
    "`beta_bm` must be an object returned from `distributional::dist_multivariate_normal`" 
    = all( 'distribution' %in% class(bm$beta_bm), family(bm$beta_bm) == 'mvnorm' )
  )
  stopifnot("`beta_bm` must refer to a single distribution" = length(bm$`beta_bm`) == 1 )
  beta_mean = mean(bm$beta_bm)[1, ]
  beta_cov  = covariance(bm$beta_bm)[[1]]
  stopifnot("Mean of `beta_unexch_bm` is not of class `numeric`" = 'numeric' %in% class(beta_mean))
  stopifnot("Covariance of `beta_unexch_bm` is not of class `matrix`" = 'matrix' %in% class(beta_cov))
  stopifnot("Mean of `beta_bm` does not match dimensions of design matrix" = length(beta_mean) == p)
  stopifnot("Mean of `beta_bm` contains missing values" = all(!is.na(beta_mean)))
  stopifnot("Covariance of `beta_bm` does not match dimensions of design matrix" = all(dim(beta_cov) == p) )
  stopifnot("Covariance of `beta_bm` contains missing values" = all(!is.na(beta_cov)))
  stopifnot("Covariance of `beta_bm` is not numerically positive definite" = all( eigen(beta_cov)$values > 1e-10 ) )
  
  stopifnot("Initial values for `beta` have not been specified" = !is.null(inits$beta))
  stopifnot("Initial values for `beta` must be a `matrix`" = 'matrix' %in% class(inits$beta) )
  stopifnot("Initial values for `beta` does not match dimensions of design matrix" = nrow(inits$beta) == p)
  stopifnot("Initial values for `beta` must be a p x K matrix" = all(nrow(inits$beta) == p, ncol(inits$beta) == K ) )
  
  hyp$beta_mean = beta_mean
  hyp$beta_cov  = beta_cov
  hyp$beta_prec = solve(beta_cov)
  
  
  
  
  ## Check precision parameter for current data: base measure and initial values
  stopifnot("`tau_bm` is not specified" = !is.null(bm$tau_bm))
  stopifnot(
    "`tau_bm` must be an object returned from `distributional::dist_gamma`" 
    = all( 'distribution' %in% class(bm$tau_bm), family(bm$tau_bm) == 'gamma' )
  )
  stopifnot("`tau_bm` must refer to a single distribution" = length(bm$tau_bm) == 1 )
  stopifnot("`tau_bm` must be an object returned from `distributional::dist_gamma`" = family(bm$tau_bm) == 'gamma' )
  
  stopifnot("Initial values for `tau have not been specified" = !is.null(inits$tau))
  stopifnot("Initial values for `tau` must be a `numeric` vector" = 'numeric' %in% class(inits$tau))
  stopifnot("Initial values for `tau` must be a K-dimensional vector" = length(inits$tau) == K)
  
  tau.params = parameters(bm$tau_bm)
  hyp$tau_shape = tau.params$shape
  hyp$tau_rate  = tau.params$rate
  
  
  ## Check concentration parameter for current data: base measure and initial values
  stopifnot("`alpha_prior` has not been specified" = !is.null(bm$alpha_prior))
  stopifnot("`alpha_prior` must refer to a single distribution" = length(bm$alpha_prior) == 1 )
  stopifnot(
    "`alpha_prior` must be an object of class `distribution` from the `distributional` package"
    = 'distribution' %in% class(bm$alpha_prior)
  )
  if ( family(bm$alpha_prior) == 'gamma' ) {
    parms = parameters(bm$alpha_prior)
    hyp$alpha_lower = 0
    hyp$alpha_upper = Inf
    hyp$alpha_shape = parms$shape
    hyp$alpha_rate  = parms$rate
  } else if ( family(bm$alpha_prior) == 'truncated' ) {
    alpha_prior_params = parameters(bm$alpha_prior)
    stopifnot( "If `alpha_prior` is truncated, its distribution must be of type `distributional::dist_gamma`" = family( parameters(bm$alpha_prior)$dist ) == 'gamma' )
    hyp$alpha_lower = alpha_prior_params$lower
    hyp$alpha_upper = alpha_prior_params$upper
    gamma_params    = parameters( parameters(bm$alpha_prior)$dist )
    hyp$alpha_shape = gamma_params$shape
    hyp$alpha_rate  = gamma_params$rate
  } else {
    stop( "`alpha_prior` must be an object from `distributional::dist_gamma` or `distributional::dist_truncated` with a `dist_gamma` distribution")
  }
  stopifnot("Initial value for `alpha` does not fit in the range of its prior" = all( inits$alpha > hyp$alpha_lower, inits$alpha < hyp$alpha_upper ) )
  
  ## Check initial values for stick-breaking variables for current data
  stopifnot("Initial values for stick-breaking weights (`v`) were not specified" = !is.null(inits$v))
  stopifnot("Initial values for stick-breaking weights (`v`) must be a `numeric` vector with elements between 0 and 1" = all( 'numeric' %in% class(inits$v), inits$v > 0, inits$v < 1 ) )
  stopifnot("Initial values for stick-breaking weights (`v`) should be of length K-1" = length(inits$v) == K-1)
  
  ## Check initial values for class assignments for current data
  stopifnot("Initial values for class assignments (`z`) were not specified" = !is.null(inits$z))
  stopifnot(
    "Initial values for class assignments (`z`) must be an n-dimensional vector of integers with elements between 1 and K (inclusive)" 
    = with(inits, all(z >= 1, z <= K, z == floor(z), z == ceiling(z), z == round(z)), !is.na(z), length(z) == n)
  )
  
  ## If external data is specified, check the appropriate base measure and initial values
  if ( !is.null(cppExtData) ) {
    
    p_unexch = ncol(cppExtData$Xunexch)
    n0       = nrow(cppExtData$Xunexch)
    
    ## Check regression coefficients for external data: base measure and initial values
    stopifnot("`beta_unexch_bm` is not specified" = !is.null(bm$beta_unexch_bm))
    stopifnot(
      "`beta_unexch_bm` must be an object returned from `distributional::dist_multivariate_normal`" 
      = all( 'distribution' %in% class(bm$beta_unexch_bm), family(bm$beta_unexch_bm) == 'mvnorm' )
    )
    stopifnot("`beta_unexch_bm` must refer to a single distribution" = length(bm$`beta_unexch_bm`) == 1 )
    beta_unexch_mean = mean(bm$beta_unexch_bm)[1, ]
    beta_unexch_cov  = covariance(bm$beta_unexch_bm)[[1]]
    stopifnot("Mean of `beta_unexch_bm` is not of class `numeric`" = 'numeric' %in% class(beta_unexch_mean))
    stopifnot("Covariance of `beta_unexch_bm` is not of class `matrix`" = 'matrix' %in% class(beta_unexch_cov))
    stopifnot("Mean of `beta_unexch_bm` does not match dimensions of design matrix; this may be because parameters in `formula` are not identifiable from the `external_data` (e.g., borrowing information from a single arm in a trial). Either use the default base measure or specify a base measure based on identifiable parameters (e.g., intercept-only model in the previous example)" = length(beta_unexch_mean) == p_unexch)
    stopifnot("Mean of `beta_unexch_bm` contains missing values" = all(!is.na(beta_unexch_mean)))
    stopifnot("Covariance of `beta_unexch_bm` does not match dimensions of design matrix" = all(dim(beta_unexch_cov) == p_unexch) )
    stopifnot("Covariance of `beta_unexch_bm` contains missing values" = all(!is.na(beta_unexch_cov)))
    stopifnot("Covariance of `beta_unexch_bm` is not numerically positive definite" = all( eigen(beta_unexch_cov)$values > 1e-10 ) )
    
    stopifnot("Initial values for `beta_unexch` have not been specified" = !is.null(inits$beta_unexch))
    stopifnot("Initial values for `beta_unexch` must be a `matrix`" = 'matrix' %in% class(inits$beta_unexch) )
    stopifnot("Initial values for `beta_unexch` does not match dimensions of design matrix" = nrow(inits$beta_unexch) == p_unexch)
    stopifnot("Initial values for `beta_unexch` must be a p_unexch x K_unexch matrix" = all(nrow(inits$beta_unexch) == p_unexch, ncol(inits$beta_unexch) == K_unexch ) )
    
    hyp$beta_unexch_mean = beta_unexch_mean
    hyp$beta_unexch_cov  = beta_unexch_cov
    hyp$beta_unexch_prec = solve( beta_unexch_cov )
    
    
    
    
    ## Check precision parameter for external data: base measure and initial values
    stopifnot("`tau_unexch_bm` is not specified" = !is.null(bm$tau_unexch_bm))
    stopifnot(
      "`tau_unexch_bm` must be an object returned from `distributional::dist_gamma`" 
      = all( 'distribution' %in% class(bm$tau_unexch_bm), family(bm$tau_unexch_bm) == 'gamma' )
    )
    stopifnot("`tau_unexch_bm` must refer to a single distribution" = length(bm$tau_unexch_bm) == 1 )
    stopifnot("`tau_unexch_bm` must be an object returned from `distributional::dist_gamma`" = family(bm$tau_unexch_bm) == 'gamma' )
    
    stopifnot("Initial values for `tau_unexch have not been specified" = !is.null(inits$tau_unexch))
    stopifnot("Initial values for `tau_unexch` must be a `numeric` vector" = 'numeric' %in% class(inits$tau_unexch))
    stopifnot("Initial values for `tau_unexch` must be a K_unexch-dimensional vector" = length(inits$tau_unexch) == K_unexch)
    
    tau_unexch_params = parameters(bm$tau_unexch_bm)
    hyp$tau_unexch_shape = tau_unexch_params$shape
    hyp$tau_unexch_rate  = tau_unexch_params$rate
    
    ## Check concentration parameter for external data: base measure and initial values
    stopifnot("`alpha_unexch_prior` has not been specified" = !is.null(bm$alpha_unexch_prior))
    stopifnot("`alpha_unexch_prior` must refer to a single distribution" = length(bm$alpha_unexch_prior) == 1 )
    stopifnot(
      "`alpha_unexch_prior` must be an object of class `distribution` from the `distributional` package"
      = 'distribution' %in% class(bm$alpha_unexch_prior)
    )
    if ( family(bm$alpha_unexch_prior) == 'gamma' ) {
      parms = parameters(bm$alpha_unexch_prior)
      hyp$alpha_unexch_lower = 0
      hyp$alpha_unexch_upper = Inf
      hyp$alpha_unexch_shape = parms$shape
      hyp$alpha_unexch_rate  = parms$rate
    } else if ( family(bm$alpha_unexch_prior) == 'truncated' ) {
      alpha_unexch_prior_params = parameters(bm$alpha_unexch_prior)
      stopifnot( "If `alpha_unexch_prior` is truncated, its distribution must be of type `distributional::dist_gamma`" = family( parameters(bm$alpha_unexch_prior)$dist ) == 'gamma' )
      hyp$alpha_unexch_lower = alpha_unexch_prior_params$lower
      hyp$alpha_unexch_upper = alpha_unexch_prior_params$lower
      gamma_params    = parameters( parameters(bm$alpha_unexch_prior)$dist )
      hyp$alpha_unexch_shape = gamma_params$shape
      hyp$alpha_unexch_rate  = gamma_params$rate
    } else {
      stop( "`alpha_unexch_prior` must be an object from `distributional::dist_gamma` or `distributional::dist_truncated` with a `dist_gamma` distribution")
    }
    stopifnot("Initial value for `alpha_unexch` does not fit in the range of its prior" = all( inits$alpha_unexch > hyp$alpha_unexch_lower, inits$alpha_unexch < hyp$alpha_unexch_upper ) )
    
    
    
    ## Check pexch prior and initial value
    stopifnot('`pexch_prior` has not been specified' = !is.null(bm$pexch_prior))
    stopifnot('`pexch_prior` must be an object of type `distribution`' = 'distribution' %in% class(bm$pexch_prior))
    params = parameters(bm$pexch_prior)
    if ( family(bm$pexch_prior) == 'beta' ) {
      hyp$pexch_shape1    = params$shape1
      hyp$pexch_shape2    = params$shape2
      hyp$pexch_lower     = 0
      hyp$pexch_upper     = 1
      hyp$pexch_priorprob = 0
    } else if ( family(bm$pexch_prior) == 'truncated' ) {
      stopifnot( "If `pexch_prior` is a truncated prior, it must be a truncated beta prior" = (family(params$dist) == 'beta'))
      beta_parameters     = parameters(params$dist)
      hyp$pexch_shape1    = beta_parameters$shape1
      hyp$pexch_shape2    = beta_parameters$shape2
      hyp$pexch_lower     = params$lower
      hyp$pexch_upper     = params$upper
      hyp$pexch_priorprob = 0
    } else if ( family(bm$pexch_prior) == 'mixture' ) {
        dists    = params$dist[[1]]
        wts      = params$w[[1]]
        famnames = sapply(dists, family)
        stopifnot(
          'If `pexch_prior` is of type `mixture`, it must be of two components: one component must be degenerate and the other must be a (possibly truncated) beta'
          = all( length(famnames) == 2, 'degenerate' %in% famnames, any( c('truncated', 'beta') %in% famnames ) )
        )
        degenindx     = which(famnames == 'degenerate')
        pointmass     = parameters(dists[[degenindx]])$x
        contdist      = dists[-degenindx][[1]]
        pointmassprob = wts[degenindx]
        
        ## Check if continuous distribution is beta
        if ( family(contdist) == 'beta' ) {
          stopifnot('If `pexch_prior` is of type `mixture`, the degenerate value must be the maximum support of the continuous distribution' = pointmass == 1)
          beta_params         = parameters(contdist)
          hyp$pexch_shape1    = beta_params$shape1
          hyp$pexch_shape2    = beta_params$shape2
          hyp$pexch_lower     = 0
          hyp$pexch_upper     = 1
          hyp$pexch_priorprob = pointmassprob
        } else {
          ## family must be truncated
          truncated_params = parameters(contdist)
          truncated_dist   = truncated_params$dist
          stopifnot( 'If `pexch_prior` is of type `mixture` and the non-degenerate component is truncated, it must be a beta distribution' = family(truncated_dist) == 'beta' )
          stopifnot( 'If `pexch_prior` is of type `mixture`, the degenerate value must be the maximum support of the continuous distribution' = pointmass == truncated_params$upper )
          beta_params = parameters(truncated_dist)
          hyp$pexch_shape1    = beta_params$shape1
          hyp$pexch_shape2    = beta_params$shape2
          hyp$pexch_lower     = max(0, truncated_params$lower)
          hyp$pexch_upper     = min(1, truncated_params$upper)
          hyp$pexch_priorprob = pointmassprob
        }
    } else if ( family(bm$pexch_prior) == 'degenerate' ) {
        hyp$pexch_shape1    = 1
        hyp$pexch_shape2    = 1
        hyp$pexch_lower     = 0
        hyp$pexch_upper     = as.numeric( parameters(bm$pexch_prior)[[1]] )
        hyp$pexch_priorprob = 1
    } else {
      stop("Prior specified on `pexch` is not supported.")
    }
    stopifnot("Initial value for `pexch` is not in prior support" = density(bm$pexch_prior, inits$pexch) > 0 )
    
    ## Check eps0 initial values
    with(inits, stopifnot( 
      '`eps0` must be a vector of length `nrow(external_data)` containing values equal to either 0 or 1' = all(eps0 == 0 | eps0 == 1, length(eps0) == n0)) 
    )
    
    ## Check initial values for stick-breaking variables for external data
    stopifnot("Initial values for stick-breaking weights (`v_unexch`) were not specified" = !is.null(inits$v_unexch))
    stopifnot("Initial values for stick-breaking weights (`v_unexch`) must be a `numeric` vector with elements between 0 and 1" = all( 'numeric' %in% class(inits$v_unexch), inits$v_unexch > 0, inits$v_unexch < 1 ) )
    stopifnot("Initial values for stick-breaking weights (`v_unexch`) should be of length K_unexch-1" = length(inits$v_unexch) == K_unexch-1)
    
    ## Check initial values for class assignments for external data
    stopifnot("Initial values for class assignments (`z0`) were not specified" = !is.null(inits$z0))
    z0    = inits$z0
    eps0  = inits$eps0
    z00   = z0[eps0==0]
    z01   = z0[eps0==1]
    stopifnot('`z0` must be a vector of length `nrow(external_data)`' = length(z0) == n0)
    stopifnot('`z0[i]` must be a value in 1:K if `eps0[i]` == 1' = all(z01 %in% 1:K))
    stopifnot('`z0[i]` must be a value in (`K`+1):(`K`+`K_unexch`) if `eps0[i]` == 0' = all(z00 %in% (K+1):(K+K_unexch)))
  }
  
  ## Return hyperparameters
  hyp
}

#' Bayesian nonparametric dynamic borrowing for time-to-event regression
#'
#' Samples from the posterior distribution of a Dirichlet process mixture model (DPMM) for time-to-event regression, using a stick-breaking representation. 
#' When `external_data` is provided, models the external data using a two-part mixture with exchangeable and nonexchangeable components, enabling dynamic borrowing.
#' Arguments given as `NULL` are optional and will have default behavior specified as below.
#'
#' @param formula a  \code{\link[stats]{formula}} giving how the survival time is related to covariates for the exchangeable data
#' @param data a \code{\link[base]{data.frame}} for the current data
#' @param formula_unexch a  \code{\link[stats]{formula}} giving how the survival time is related to covariates for the unexchangeable data; this may be different form `formula` when, e.g., borrowing information from a control arm only.
#' @param external_data optional \code{\link[base]{data.frame}} for the external data. If `NULL`, a standard DPMM is fit to `data`.
#' @param K integer with value \eqn{\ge 2} giving truncation point of the stick-breaking representation for the exchangeable population
#' @param prior_ess a value giving the prior effective sample size for the base measure; ignored if the base measure is specified by the user
#' @param K_unexch truncation point of the stick-breaking representation for the nonexchangeable population; required to be an integer \eqn{\ge 2} if `external_data` is not `NULL`
#' @param bm an optional named `list` containing base measures and priors for parameters with the following elements:
#' \itemize{
#'  \item `beta_bm`: either `NULL` or a \code{\link[distributional]{dist_multivariate_normal}} base measure for regression coefficients of the exchangeable population; if `NULL`, defaults to a multivariate normal prior with hyperparameters based on the maximum likelihood estimate and inverse information matrix (rescaled to `prior_ess` observations) using `data`
#'  \item `tau_bm`: either `NULL` or a \code{\link[distributional]{dist_gamma}} base measure on the precision parameter for the exchangeable population; if `NULL`, defaults to a gamma prior with shape parameter \eqn{\ge} 1.1 that most closely approximates the log-normal distribution implied by the maximum likelihood analysis with asymptotic variance rescaled to `nsubj` observations using `data`
#'  \item `alpha_prior`: either `NULL` or a \code{\link[distributional]{dist_gamma}} prior or a truncated gamma prior via \code{\link[distributional]{dist_truncated}} for the concentration parameter of the exchangeable DPMM. if `NULL`, defaults to a gamma prior with shape parameter 8 and rate parameter 2
#'  \item `beta_unexch_bm`: either `NULL` or a \code{\link[distributional]{dist_multivariate_normal}} base measure for regression coefficients of the nonexchangeable population; required if `external_data` is not `NULL`; if `NULL`, defaults to a multivariate normal prior with hyperparameters based on the maximum likelihood estimate and inverse information matrix (rescaled to `nsubj` observations) using `external_data`
#'  \item `tau_unexch_bm`: either `NULL` or a \code{\link[distributional]{dist_gamma}} base measure for the precision of the nonexchangeable population; required if `external_data` is not `NULL`; if `NULL`, defaults to a gamma prior with shape parameter \eqn{\ge} 1.1 that most closely approximates the log-normal distribution implied by the maximum likelihood analysis with asymptotic variance rescaled to `nsubj` observations using `external_data`
#'  \item `alpha_unexch_prior`: either `NULL` or a \code{\link[distributional]{dist_gamma}} prior or a truncated gamma prior via \code{\link[distributional]{dist_truncated}} for the concentration parameter of the nonexchangeable DPMM; if `NULL`, defaults to a gamma prior with shape parameter 8 and rate parameter 2
#'  \item `pexch_prior`: either `NULL`, a \code{\link[distributional]{dist_beta}}, a truncated beta prior via \code{\link[distributional]{dist_truncated}}, a fixed value via \code{\link[distributional]{dist_degenerate}}, or a mixture of a degenerate distribution and a beta or truncated beta prior via \code{\link[distributional]{dist_mixture}} prior for the probability that an external observation is exchangeable; if `NULL`, defaults a 0.5-0.5 mixture of a uniform (Beta(1,1)) prior and a point mass at 1. If the prior is a mixture, the degenerate distribution must line up with the maximum possible value of `pexch`.
#' }
#' @param nburnin number of MCMC burn-in samples
#' @param nsamples number of desired posterior samples after burn-in (and thinning)
#' @param thin thinning parameter for MCMC sampling
#' @param inits an optional `list` containing starting values for any of the following:
#' \itemize{
#'   \item `beta`: \eqn{K \times p} `matrix` giving cluster regression coefficients for the exchangeable population; defaults to random draws from a multivariate normal centered at the MLE using `data` with covariance equal to the inverse Fisher Information
#'   \item `tau`: \eqn{K}-`vector` giving cluster precision parameters for the exchangeable population; defaults to random draws from a log-normal centered at the MLE using `data` with covariance equal to the inverse Fisher Information
#'   \item `alpha`: Dirichlet process concentration parameter for exchangeable population; set to the prior mean implied by `alpha_prior`
#'   \item \code{v}: \eqn{(K - 1)}-dimensional vector giving stick-breaking variables for the exchangeable population; defaults to \eqn{(K - 1)} samples from a Beta(1, \code{alpha}) distribution
#'   \item `pexch`: probability that a subject in `external_data` is exchangeable; defaults to prior mean of `pexch_prior`; ignored if `external_data` is `NULL`
#'   \item `eps0`: \eqn{n_0 \times K_{unexch}} matrix of exchangeability indicators for external data; default value is based on highest likelihood from initial values; ignored if `external_data` is `NULL`
#'   \item `beta_unexch`: \eqn{K_{unexch} \times p} `matrix` giving cluster regression coefficients for unexchangeable population; default behavior is analogous to `beta` but using `external_data`; ignored if `external_data` is `NULL`
#'   \item `tau_unexch`: \eqn{K_{unexch}}-`vector` giving cluster precision parmaeters for unexchangeable population; default behavior is analogous to `tau` but using `external_data`; ignored if `external_data` is `NULL`
#'   \item `alpha_unexch`: mean of `alpha_unexch_prior`; ignored if `external_data` is `NULL`
#'   \item `v_unexch`: \eqn{(K_{unexch} - 1)}-vector giving stick-breaking variables for unexchangeable population; defaults to samples from Beta(1, alpha_unexch); ignored if `external_data` is `NULL`
#'   \item `z0`: a \eqn{n_0}-`vector` consisting of starting values for latent clusters for the external data; default is based on likelihood using initial parameters; ignored if `external_data` is `NULL`; if `eps0 == 1` then `z0` should be between `1:K` otherwise between `1:K_unexch`. Defaults to the maximum log likelihood among the `K` clusters using `beta_init`, `tau_init`, `beta_unexch_init`, and `tau_unexch_init`. Note: if `z0` is specified then `eps0` must also be specified.
#' }
#'
#' @return A `list` (class `"bnpdb_tte_conditional"`) with the following components:
#' \itemize{
#'   \item `beta`: a \eqn{p \times K \times M} `array` of regression coefficients for the exchangeable population
#'   \item `sigma`: an \eqn{M \times K} `matrix` of standard deviations for the exchangeable population
#'   \item `alpha`: an `M`-dimensional `vector` of concentration parameters for the exchangeable DPMM
#'   \item `w`: a \eqn{M \times K} `matrix` of stick-breaking weights for the exchangeable population
#'   \item `pexch`: an `M`-dimensional `vector` of posterior draws of exchangeability probability; only present if `external_data` is not `NULL`
#'   \item `eps0`: a \eqn{M \times n_0} `matrix` of exchangeability indicators for external data; only present if `external_data` is not `NULL`
#' }
#' 
#' @examples
#' ## Set seed for reproducibility
#' set.seed(741)
#' 
#' ## Generate current and historical data sets
#' n     <- 100
#' n0    <- 100
#' N     <- n + n0
#' pexch <- 0.80
#' beta  <- cbind(c(1, 1, 0), c(-1, -1, 0))
#' sigma <- c(0.5, 0.75)
#' a     <- rbinom(N, 1, 0.5)      ## tretment indicator
#' x     <- rbinom(N, 1, 0.5)      ## binary covariate
#' eps0  <- rbinom(n0, 1, pexch)   ## exchangeability indicator
#' eps   <- c(rep(1, n), 2 - eps0) ## 1 = exch; 2 = unexch
#' X     <- cbind(1, a, x)
#' Mu    <- unlist( 
#'   lapply(1:N, function(i) { as.numeric(X[i, ] %*% beta[, eps[i]] ) } ) 
#' )
#' logt     <- rnorm(N, Mu, sigma[eps])
#' logc     <- rnorm(N, mean(Mu) + 0.50, sd = 0.25)
#' event    <- logt <= logc
#' logy     <- ifelse(event, logt, logc)
#' dat      <- data.frame(y = exp(logy), event = event, a = a, x = x)
#' curdata  <- dat[1:n, ]
#' histdata <- dat[-(1:n), ]
#' 
#' ## Times for marginalization
#' tau   = max(curdata$y[curdata$event == 1])
#' times = seq(0.01, tau, length.out = 5)
#' 
#' ## Fit intercept only (no borrowing)
#' fit.intonly = tte_reg_bnp(
#'   survival::Surv(y, event) ~ 1, data = curdata, K = 5
#'   , nburnin = 0, nsamples = 100
#' )
#' 
#' ## Fit treatment only (with borrowing)
#' fit.strataonly = tte_reg_bnp(
#'   survival::Surv(y, event) ~ a, data = curdata, K = 5
#'   , external_data = histdata, K_unexch = 5
#'   , nburnin = 0, nsamples = 100
#' )
#' 
#' 
#' ## Fit treatment and binary covariate (with borrowing)
#' fit.anova_ddp = tte_reg_bnp(
#'   survival::Surv(y, event) ~ a + x, data = curdata, K = 5
#'   , external_data = histdata, K_unexch = 5
#'   , nburnin = 0, nsamples = 100
#' )
#' 
#' ## Fit models stratified by arm (no borrowing)--results in a list.
#' fit.stratified = lapply(
#'   0:1, function(arm) {
#'     tte_reg_bnp(
#'       survival::Surv(y, event) ~ x
#'       , data = curdata[curdata$a == arm, ]
#'       , K = 5, K_unexch = 5
#'       , nburnin = 0, nsamples = 100
#'     )
#'   }
#' )
#' 
#' ## Summarize conditional estimates--only reports identifiable quantities
#' summary(fit.strataonly, mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)))
#' summary(fit.anova_ddp)
#' lapply(fit.stratified, summary)
#'
#' @export
tte_reg_bnp = function(
    formula, data, K
    , formula_unexch = formula, external_data = NULL, K_unexch = NULL
    , prior_ess = 5
    , bm = list()
    , inits = list()
    , nburnin = 2000, nsamples = 10000, thin = 1
) {
  stopifnot("`formula` must be an object of type formula" = inherits(formula, "formula"))
  stopifnot("`data` must be an object of type data.frame" = is.data.frame(data))
  stopifnot("`K` must be an integer >= 2`" = is.numeric(K), K %% 1 == 0, K >= 2)
  
  ## Fit current data MLE; extract data components
  fit     = survival::survreg(formula, data = data, dist = 'lognormal')
  mle     = c(coef(fit), 'Log(scale)' = log(fit$scale))
  if ( any(is.na(mle) ) )
    stop("`formula` contains unidentifiable parameters based on `data`. Try reparamterizing the model so all parameters are identifiable")
  vcov       = vcov(fit)
  X          = model.matrix(fit, data = data)
  y          = fit$y[, 1]
  event      = fit$y[, 2]
  n          = length(y)
  p          = ncol(X)
  mle0       = NULL
  vcov0      = NULL
  n0         = NULL
  cppData    = list(y = y, event = event, X = X)
  cppExtData = NULL
  
  ## Get data for C++ code
  if ( !is.null(external_data) ) {
    stopifnot("`formula_unexch` must be an object of type formula" = inherits(formula_unexch, "formula"))
    stopifnot("`external_data` must be an object of type data.frame" = is.data.frame(external_data))
    stopifnot("`K_unexch` must be an integer >= 2`" = is.numeric(K_unexch), K_unexch %% 1 == 0, K_unexch >= 2)
    fit0  = survival::survreg(formula_unexch, data = external_data, dist = 'lognormal')
    mle0  = c(coef(fit0), 'Log(scale)' = log(fit$scale))
    if ( any(is.na(mle0) ) )
      stop("`formula_unexch` contains unidentifiable parameters based on `external_data`. This can happen, e.g., when borrowing information from only a control arm. Try reparamterizing the model via `formula_unexch` so all parameters are identifiable")
    vcov0 = vcov(fit0)
    
    cppExtData$Xunexch = model.matrix(fit0, data = external_data)
    cppExtData$y       = log( fit0$y[, 1] )
    cppExtData$event   = fit0$y[, 2]
    cppExtData$X       = model.matrix(fit, data = external_data)
    n0 = length(cppExtData$y)
    p0 = ncol(cppExtData$Xunexch)
  }
  
  ## Augment base measure with defaults if components are unspecified
  bm = get_bm_tte(bm, n, mle, vcov, mle0, vcov0, n0, prior_ess)
  
  ## Obtain initial values (if user did not specify) 
  inits = tte_inits(inits, bm, cppData, K, cppExtData, K_unexch)
  
  ## Check base measure and initial values; return hyperparameters
  hyperparams = tte_check_all_and_export_hyperparams(cppData, K, bm, inits, cppExtData, K_unexch)
  
  ## Put cpp data time on log scale for sampling
  cppData$y = log(cppData$y)
  
  
  ## If external data is not specified, proceed to C++ sampling of DPMM using
  ## only specified `data`
  
  if (is.null(external_data)) {
    smpl = tte_reg_dpmm_cpp(cppData, hyperparams, inits, nsamples, nburnin, thin)
    
  } else {
    ## If `external_data` is specified, do the same checks and initial values
    ## as above before proceeding to MCMC sampling
    
    ## Obtain posterior sampling using BNPDB
    smpl = tte_reg_bnpdb_cpp(cppData, cppExtData, hyperparams, inits, nsamples, nburnin, thin)
    
    colnames(smpl$eps0) = paste0('eps0[', 1:nrow(cppExtData$X), ']')
    colnames(smpl$pexch) = "pexch"
  }
  
  rownames(smpl$beta) = colnames(X)
  colnames(smpl$beta) = paste0("K=", 1:K)
  colnames(smpl$sigma) = paste0("sigma[", 1:K, "]")
  colnames(smpl$w) = paste0("w[", 1:K, "]")
  colnames(smpl$alpha) = "alpha"
  
  attr(smpl, "formula") = formula
  attr(smpl, "inits") = inits
  attr(smpl, "data") = data
  attr(smpl, "external_data") = external_data
  attr(smpl, "basemeasure") = bm
  attr(smpl, "datatype") = "tte"
  attr(smpl, "type") = "conditional parameters"
  class(smpl) = c("bnpdb_tte_conditional", class(smpl))
  
  smpl
}
