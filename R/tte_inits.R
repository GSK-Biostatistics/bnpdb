


##--------------------------------------------------------------------------
## This file contains functions to check the initial values for 
## continuous DPMM models, using sensible defaults when the user has left
## them unspecified.
##
## This file contains helper functions, and is not part of the main output
## of bnpdb.
##--------------------------------------------------------------------------


#' Get and check initial values for internal data / exchangeable parameters (if unspecified)
#' @param inits named list of initial values
#' @param fit fitted survival regression w/ log-normal distribution to internal data
#' @param K number of components in mixture
#' @param alpha_prior prior on concentration parameter
#' @return named list giving initial values pertaining to internal / exchangeable parameters
#' 
#' @importFrom mvtnorm rmvnorm
#' 
#' @keywords internal
#' @noRd
tte_get_inits_internal = function(formula, data, fit, K, alpha_prior, inits = list()) {
  # Construct design matrix and survival outcome directly
  X     = model.matrix(formula, data)
  y     = fit$y[, 1]
  event = fit$y[, 2]
  n     = nrow(X)
  p     = ncol(X)
  
  # Extract MLE and covariance matrix
  mle  = c(coef(fit), 'Log(scale)' = log(fit$scale))
  vcov = vcov(fit)
  
  # Use user-supplied inits or generate sensible defaults
  user_inits = names(inits)
  
  if (!('beta' %in% user_inits)) {
    inits$beta = t(mvtnorm::rmvnorm(K, mle[1:p], vcov[1:p, 1:p, drop = FALSE]))
  }
  
  if (!('tau' %in% user_inits)) {
    inits$tau = rlnorm(K, meanlog = mle[p + 1], sdlog = sqrt(vcov[p + 1, p + 1]))
  }
  
  if (!('z' %in% user_inits)) {
    Mu    = X %*% inits$beta
    sigma = 1 / sqrt(inits$tau)
    logf  = matrix(nrow = n, ncol = K)
    for (k in 1:K) {
      logf[, k] = ifelse(
        event,
        dnorm(y, Mu[, k], sigma[k], log = TRUE),
        pnorm(y, Mu[, k], sigma[k], lower.tail = FALSE, log.p = TRUE)
      )
    }
    inits$z = apply(logf, 1, which.max)
  }
  
  if (!('alpha' %in% user_inits)) {
    inits$alpha = mean(alpha_prior)
  }
  
  if (!('v' %in% user_inits)) {
    inits$v = rbeta(K - 1, 1, inits$alpha)
  }
  
  # Validate inits
  stopifnot(
    "Initial value for beta must be a p x K `matrix` with real values" = {
      is.matrix(inits$beta) &&
        nrow(inits$beta) == p &&
        ncol(inits$beta) == K &&
        all(is.finite(inits$beta))
    },
    "Initial value for tau must be a `numeric` vector of length `K` with only positive elements" = {
      is.numeric(inits$tau) &&
        length(inits$tau) == K &&
        all(is.finite(inits$tau) & inits$tau > 0)
    },
    "Initial value for z must be a vector of length `nrow(data)` with integer values 1:K" = {
      is.numeric(inits$z) &&
        length(inits$z) == n &&
        all(inits$z %% 1 == 0 & inits$z >= 1 & inits$z <= K)
    },
    "Initial value for alpha must be a positive scalar obeying the bounds of its (possibly truncated) gamma prior" = {
      length(inits$alpha) == 1 &&
        inits$alpha >= unlist(support(alpha_prior))["lim1"] &&
        inits$alpha <= unlist(support(alpha_prior))["lim2"] &&
        is.finite(inits$alpha)
    },
    "Initial value for v must be a vector of length `K-1` consisting of values between 0 and 1 exclusive" = {
      length(inits$v) == K - 1 &&
        all(inits$v > 0 & inits$v < 1)
    }
  )
  
  return(inits)
}

#' Get initial values for external data / unexchangeable parameters (if unspecified)
#' @param inits named list of initial values
#' @param fit fitted survival regression w/ log-normal distribution to internal data
#' @param K number of components in mixture
#' @param alpha_unexch_prior prior on concentration parameter
#' @return named list giving initial values pertaining to internal / exchangeable parameters
#' @keywords internal
#' @noRd
tte_get_inits_external = function(formula, external_data, fit0, K_unexch, alpha_unexch_prior, inits = list()) {
  K      = ncol(inits$beta)
  X0     = model.matrix(fit0, external_data)
  n0     = nrow(X0)
  p      = ncol(X0)
  mle0   = c(coef(fit0), 'Log(scale)' = log(fit0$scale))
  vcov0  = vcov(fit0)
  event0 = fit0$y[, 2]
  y0     = fit0$y[, 1]
  
  ##---------------------------------------------------------------------------
  ## If user did not specify initial values, use defaults
  ##---------------------------------------------------------------------------
  user_inits = names(inits)
  if ( !( 'beta_unexch' %in% user_inits ))
    inits$beta_unexch = t( mvtnorm::rmvnorm(K_unexch, mle0[1:p], vcov0[1:p, 1:p, drop = FALSE]) )
  if ( !( 'tau_unexch' %in% user_inits ))
    inits$tau_unexch = rlnorm(K_unexch, meanlog = mle0[p+1], sdlog = sqrt(vcov0[p+1,p+1]))
  if ( !('eps0' %in% user_inits) & ('z0' %in% user_inits) ) {
    stop("Initial values were specified for z0 but not eps0. Either specify an initial value for eps0 as well as z0 or leave both unspecified.")
  }
  sigma  = inits$tau
  sigma0 = inits$tau_unexch
  if ( ('eps0' %in% user_inits) & !('z0' %in% user_inits) ) {
    message('Initial values were specified for eps0 but not z0. Initial values for z0 will be set based on log likelihood contributions and iniital values for beta, tau, beta_unexch, and tau_unexch')
    Mu      = X0 %*% inits$beta
    Mu0     = X0 %*% inits$beta_unexch
    sigma0  = 1 / sqrt(inits$tau_unexch)
    logf    = matrix(nrow = n0, ncol = K)
    logf0   = matrix(nrow = n0, ncol = K_unexch)
    for ( k in 1:K )
      logf[, k] = ifelse(event0, dnorm(y0, Mu[, k], sigma[k], log = TRUE), pnorm(y0, Mu[, k], sigma[k], lower.tail = FALSE, log.p = TRUE) )
    for ( k in 1:K_unexch )
      logf0[, k] = ifelse(event0, dnorm(y0, Mu0[, k], sigma0[k], log = TRUE), pnorm(y0, Mu0[, k], sigma0[k], lower.tail = FALSE, log.p = TRUE) )
    
    inits$z0 = rep(NA, n0)
    inits$z0[inits$eps0 == 1] = apply(logf[inits$eps0 == 1, ], 1, which.max)
    inits$z0[inits$eps0 == 0] = apply(logf0[inits$eps0 == 0, ], 1, which.max)
  }
  if ( !( 'eps0' %in% user_inits | 'z0' %in% user_inits ) ) {
    Mu      = X0 %*% inits$beta
    Mu0     = X0 %*% inits$beta_unexch
    sigma0  = 1 / sqrt(inits$tau_unexch)
    logf    = matrix(nrow = n0, ncol = K)
    logf0   = matrix(nrow = n0, ncol = K_unexch)
    for ( k in 1:K )
      logf[, k] = ifelse(event0, dnorm(y0, Mu[, k], sigma[k], log = TRUE), pnorm(y0, Mu[, k], sigma[k], lower.tail = FALSE, log.p = TRUE) )
    for ( k in 1:K_unexch )
      logf0[, k] = ifelse(event0, dnorm(y0, Mu0[, k], sigma0[k], log = TRUE), pnorm(y0, Mu0[, k], sigma0[k], lower.tail = FALSE, log.p = TRUE) )
    indx0        = apply(cbind(logf, logf0), 1, which.max)
    inits$eps0   = ifelse(indx0 <= K, 1, 0)
    inits$z0     = ifelse(inits$eps0 == 1, indx0, indx0 - K)
  }
  if ( !('alpha_unexch' %in% user_inits) )
    inits$alpha_unexch = mean(alpha_unexch_prior)
  
  if ( !('v_unexch' %in% user_inits) )
    inits$v_unexch = rbeta(K_unexch - 1, 1, inits$alpha_unexch)
  
  if ( !('pexch' %in% user_inits) )
    inits$pexch = mean(inits$eps0)
  
  ##---------------------------------------------------------------------------
  ## Check initial values
  ##---------------------------------------------------------------------------
  stopifnot( 
    "Initial value for beta_unexch must be a p x K_unexch `matrix` with real values" = {
      'matrix' %in% class(inits$beta_unexch) 
      nrow(inits$beta_unexch) == p
      ncol(inits$beta_unexch) == K_unexch
      !any(is.null(inits$beta_unexch), is.na(inits$beta_unexch), is.infinite(inits$beta_unexch))
    })
  stopifnot( 
    "Initial value for tau_unexch must be a `numeric` vector of length `K_unexch` with only positive elements" = {
      'numeric' %in% class(inits$tau_unexch) 
      length(inits$tau_unexch) == K_unexch
      !any(is.null(inits$tau_unexch), is.na(inits$tau_unexch), is.infinite(inits$tau_unexch), inits$tau_unexch <= 0)
    })
  stopifnot( 
    "Initial value for eps0 must be a vector of length `nrow(external_data)` with integer values 0:1 with 0 = unexchangeable and 1 = exchangeable" = { 
      any( 'numeric' %in% class(inits$eps0), 'integer' %in% class(inits$eps0) )
      length(inits$eps0) == n0
      !any(is.null(inits$eps0), is.na(inits$eps0), is.infinite(inits$eps0), inits$eps0 < 0, inits$eps0 > 1)
      all(inits$eps0 %% 1 == 0)
    })
  stopifnot( 
    "Initial value for z0 must be a vector of length `nrow(external_data)` with integer values 1:K_unexch if eps0 == 0 [unexchangeable] and 1:K if eps0 == 1 [exchangeable]" = { 
      any( 'numeric' %in% class(inits$z0), 'integer' %in% class(inits$z0) )
      length(inits$z0) == n0
      !any(is.null(inits$z0), is.na(inits$z0), is.infinite(inits$z0), inits$z0 <= 0)
      all(inits$z0 %% 1 == 0)
      inits$z0[inits$eps0 == 1] >= 1; inits$z0[inits$eps0 == 1] <= K
      inits$z0[inits$eps0 == 0] >= 1; inits$z0[inits$eps0 == 0] <= K_unexch
    })
  stopifnot( 
    "Initial value for alpha_unexch must be a positive scalar obeying the bounds of its (possibly truncated) gamma prior" = { 
      length(inits$alpha_unexch) == 1
      inits$alpha_unexch >= unlist(support(alpha_unexch_prior))['lim1']
      inits$alpha_unexch <= unlist(support(alpha_unexch_prior))['lim2']
      is.finite(inits$alpha_unexch); !is.na(inits$alpha_unexch)
    })
  stopifnot( 
    "Initial value for v_unexch must be a vector of length `K_unexch-1` consisting of values between 0 and 1 exclusive" = { 
      length(inits$v_unexch) == K_unexch-1; min(inits$v_unexch) > 0; max(inits$v_unexch) < 1
    })
  
  stopifnot(
    "Initial value for `pexch` must be a scalar between 0 and 1 (exclusive)" = {
      length(inits$pexch) == 1; inits$pexch > 0; inits$pexch < 1
    }
  )
  
  ## Return initial values
  return(inits)
}




