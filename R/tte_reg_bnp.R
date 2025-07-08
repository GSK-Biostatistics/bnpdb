#' Bayesian nonparametric dynamic borrowing for time-to-event regression
#'
#' Samples from the posterior distribution of a Dirichlet process mixture model (DPMM) for time-to-event regression, using a stick-breaking representation. 
#' When `external_data` is provided, models the external data using a two-part mixture with exchangeable and nonexchangeable components, enabling dynamic borrowing.
#' Arguments given as `NULL` are optional and will have default behavior specified as below.
#'
#' @include inputs_continuous.R
#' @include tte_inits.R
#'
#' @param formula a  \code{\link[stats]{formula}} giving how the survival time is related to covariates
#' @param data a \code{\link[base]{data.frame}} for the current data
#' @param external_data optional \code{\link[base]{data.frame}} for the external data. If `NULL`, a standard DPMM is fit to `data`.
#' @param K integer with value \eqn{\ge 2} giving truncation point of the stick-breaking representation for the exchangeable population
#' @param K_unexch truncation point of the stick-breaking representation for the nonexchangeable population; required to be an integer \eqn{\ge 2} if `external_data` is not `NULL`
#' @param beta_bm either `NULL` or a \code{\link[distributional]{dist_multivariate_normal}} base measure for regression coefficients of the exchangeable population; if `NULL`, defaults to a multivariate normal prior with hyperparameters based on the maximum likelihood estimate and inverse information matrix (rescaled to 5 events) using `data`.
#' @param tau_bm a \code{\link[distributional]{dist_gamma}} base measure on the precision parameter for the exchangeable population; if `NULL`, defaults to a gamma prior that most closely approximates the log-normal distribution implied by the maximum likelihood analysis with variance rescaled to 5 events using `data`
#' @param alpha_prior a \code{\link[distributional]{dist_gamma}} prior or a truncated gamma prior via \code{\link[distributional]{dist_truncated}} for the concentration parameter of the exchangeable DPMM
#' @param beta_unexch_bm either `NULL` or a \code{\link[distributional]{dist_multivariate_normal}} base measure for regression coefficients of the nonexchangeable population; required if `external_data` is not `NULL`; if `NULL`, defaults to a multivariate normal prior with hyperparameters based on the maximum likelihood estimate and inverse information matrix (rescaled to 5 events) using `external_data`
#' @param tau_unexch_bm a \code{\link[distributional]{dist_gamma}} base measure for the precision of the nonexchangeable population; required if `external_data` is not `NULL`; if `NULL`, defaults to a gamma prior that most closely approximates the log-normal distribution implied by the maximum likelihood analysis with variance rescaled to 5 events using `data`
#' @param alpha_unexch_prior a \code{\link[distributional]{dist_gamma}} prior or a truncated gamma prior via \code{\link[distributional]{dist_truncated}} for the concentration parameter of the nonexchangeable DPMM; required if `external_data` is not `NULL`
#' @param pexch_prior a \code{\link[distributional]{dist_beta}} prior for the probability that an external observation is exchangeable; required if `external_data` is not `NULL`
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
#'   \item `eps0`: \eqn{n_0 \times K_unexch} matrix of exchangeability indicators for external data; default value is based on highest likelihood from initial values; ignored if `external_data` is `NULL`
#'   \item `beta_unexch`: \eqn{K_unexch \times p} `matrix` giving cluster regression coefficients for unexchangeable population; default behavior is analogous to `beta` but using `external_data`; ignored if `external_data` is `NULL`
#'   \item `tau_unexch`: \eqn{K_unexch}-`vector` giving cluster precision parmaeters for unexchangeable population; default behavior is analogous to `tau` but using `external_data`; ignored if `external_data` is `NULL`
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
    , external_data = NULL, K_unexch = NULL
    , beta_bm = NULL, tau_bm = NULL
    , alpha_prior = dist_truncated(dist_gamma(12, 2), 0, Inf)
    , beta_unexch_bm = NULL, tau_unexch_bm = NULL
    , alpha_unexch_prior = dist_truncated(dist_gamma(12, 2), 0, Inf)
    , pexch_prior = dist_beta(0.5, 0.5)
    , inits = list()
    , nburnin = 2000, nsamples = 10000, thin = 1
) {
  stopifnot("`formula` must be an object of type formula" = inherits(formula, "formula"))
  stopifnot("`data` must be an object of type data.frame" = is.data.frame(data))
  stopifnot("`K` must be an integer >= 2`" = is.numeric(K), K %% 1 == 0, K >= 2)
  
  ## Fit current data MLE
  fit    = survival::survreg(formula, data = data, dist = 'lognormal')
  X      = model.matrix(fit, data = data)
  y      = fit$y[, 1]
  event  = fit$y[, 2]
  
  ## Get data for C++ code
  data_cpp = list(y = log(y), event = event, X = X)
  
  ## Check priors / base measures and obtain hyperparameters from priors as list
  priors      = check_prior_list_continuous(fit, alpha_prior, beta_bm, tau_bm, X)
  hyperparams = priors$hyperparams
  
  ## Obtain (if user did not specify) and check initial values
  inits = tte_get_inits_internal(formula, data, fit, K, alpha_prior, inits)
  
  if (is.null(external_data)) {
    
    ## If external data is not specified, proceed to C++ sampling of DPMM using
    ## only specified `data`
    smpl = tte_reg_dpmm_cpp(data_cpp, hyperparams, inits, nsamples, nburnin, thin)
    
  } else {
    
    ## If `external_data` is specified, do the same checks and initial values
    ## as above before proceeding to MCMC sampling
    stopifnot("If specified, `external_data` must be an object of type data.frame" = is.data.frame(external_data))
    stopifnot("If `external_data` is not NULL, K_unexch` must be an integer >= 2`" = is.numeric(K_unexch), K_unexch %% 1 == 0, K_unexch >= 2)
    
    ## Fit MLE to external_data
    fit0   = survival::survreg(formula, data = external_data, dist = 'lognormal')
    X0     = model.matrix(fit0, data = external_data)
    y0     = fit0$y[, 1]
    event0 = fit0$y[, 2]
    
    ## Get C++ version of external data
    external_data_cpp = list(y = log(y0), event = event0, X = X0)
    
    ## Obtain (if user did not specify) and check initial values
    inits = tte_get_inits_external(formula, external_data, fit0, K_unexch, alpha_unexch_prior, inits)
    
    ## Check priors / base measures and add to hyperparmaeters for the current data
    priors_ext  = check_prior_list_continuous(fit0, alpha_unexch_prior, beta_unexch_bm, tau_unexch_bm, X0, external = TRUE)
    hyperparams = c(hyperparams, priors_ext$hyperparams)
    stopifnot( "`pexch_prior` must be an object from `distributional::dist_beta()`" = 'beta' == family(pexch_prior))
    hyperparams$pexch_shape1 = distributional::parameters(pexch_prior)[['shape1']]
    hyperparams$pexch_shape2 = distributional::parameters(pexch_prior)[['shape2']]
    
    ## Obtain posterior sampling using BNPDB
    smpl = tte_reg_bnpdb_cpp(data_cpp, external_data_cpp, hyperparams, inits, nsamples, nburnin, thin)
    
    colnames(smpl$eps0) = paste0('eps0[', 1:nrow(X0), ']')
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
  attr(smpl, "basemeasure") = list(
    alpha_prior = priors$alpha_prior
    , beta_bm = beta_bm
    , tau_bm = tau_bm
    , alpha_unexch_prior = if (!is.null(external_data)) alpha_unexch_prior else NULL
    , beta_unexch_bm     = if (!is.null(external_data)) beta_unexch_bm else NULL
    , tau_unexch_bm      = if (!is.null(external_data)) tau_unexch_bm else NULL
    , basemeasure        = if (!is.null(external_data)) c(priors, priors_ext) else priors
  )
  attr(smpl, "datatype") = "tte"
  attr(smpl, "type") = "conditional parameters"
  class(smpl) = c("bnpdb_tte_conditional", class(smpl))
  
  smpl
}
