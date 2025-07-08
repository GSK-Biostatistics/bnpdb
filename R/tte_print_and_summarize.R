#' Summarize a bnpdb time-to-event regression model with conditional parameters
#' 
#' Summarize a bnpdb time-to-event regression model with conditional parameters
#' 
#' @param object An object of class `bnpdb_tte_conditional` obtained by calling \code{\link{tte_reg_bnp}}
#' @param ... arguments to pass onto \code{\link[posterior]{summarize_draws}}
#' 
#' @return an object of class 
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
#' ## Fit models stratified by arm (no borrowing)
#' fit.stratified = lapply(
#'   0:1, function(arm) {
#'     tte_reg_bnp(
#'       survival::Surv(y, event) ~ x
#'       , data = curdata[curdata$a == arm, ]
#'       , external_data = histdata[histdata$a == arm, ]
#'       , K = 5, K_unexch = 5
#'       , nburnin = 0, nsamples = 100
#'     )
#'   }
#' )
#' 
#' ## Summarize conditional estimates--only reports identifiable quantities
#' summary(fit.strataonly, mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)))
#' summary(fit.anova_ddp)
#'
#' @importFrom posterior summarize_draws quantile2
#' 
#' @export
#' @method summary bnpdb_tte_conditional
summary.bnpdb_tte_conditional <- function(object, ...) {
  if (is.null(attr(object, "external_data"))) {
    temp <- object[c("alpha", "w")]
  } else {
    temp <- object[c("alpha", "w", "pexch")]
  }
  temp <- do.call("cbind", temp)
  out <- posterior::summarize_draws(temp, ...)
  class(out) <- c("bnpdb_tte_conditional_summary", class(out))
  out
}



#' Summarize a bnpdb time-to-event regression model with marginal parameters
#' 
#' Summarize a bnpdb time-to-event regression model with marginal parameters
#' 
#' @param object An object of class `bnpdb_tte_marginal` obtained from calling \code{\link{tte_marginalize}} (or a `list` of such objects).
#' @param ... arguments to pass onto \code{\link[posterior]{summarize_draws}}
#' 
#' @return a named list consisting of the following:
#' \itemize{
#'   \item `beta`: \eqn{K \times p} `matrix` giving cluster regression coefficients for the exchangeable population; defaults to random draws from a multivariate normal centered at the MLE using `data` with covariance equal to the inverse Fisher Information
#'   \item `tau`: \eqn{K}-`vector` giving cluster precision parameters for the exchangeable population; defaults to random draws from a log-normal centered at the MLE using `data` with covariance equal to the inverse Fisher Information
#'   \item `alpha`: Dirichlet process concentration parameter for exchangeable population; set to the prior mean implied by `alpha_prior`
#'   \item `v`: \eqn{(K-1)} `vector` giving stick-breaking variables for the exchangeable population; defaults to \eqn{(K - 1)} samples from a Beta(1, alpha) distribution
#'   \item `pexch`: probability that a subject in `external_data` is exchangeable; defaults to prior mean of `pexch_prior`; ignored if `external_data` is `NULL`
#'   \item `eps0`: \eqn{n_0 \times K_unexch} matrix of exchangeability indicators for external data; default value is based on highest likelihood from initial values; ignored if `external_data` is `NULL`
#'   \item `beta_unexch`: \eqn{K_unexch \times p} `matrix` giving cluster regression coefficients for unexchangeable population; default behavior is analogous to `beta` but using `external_data`; ignored if `external_data` is `NULL`
#'   \item `tau_unexch`: \eqn{K_unexch}-`vector` giving cluster precision parmaeters for unexchangeable population; default behavior is analogous to `tau` but using `external_data`; ignored if `external_data` is `NULL`
#'   \item `alpha_unexch`: mean of `alpha_unexch_prior`; ignored if `external_data` is `NULL`
#'   \item `v_unexch`: \eqn{(K_{unexch} - 1)}-vector giving stick-breaking variables for unexchangeable population; defaults to samples from Beta(1, alpha_unexch); ignored if `external_data` is `NULL`
#'   \item `z0`: a \eqn{n_0}-`vector` consisting of starting values for latent clusters for the external data; default is based on likelihood using initial parameters; ignored if `external_data` is `NULL`; if `eps0 == 1` then `z0` should be between `1:K` otherwise between `1:K_unexch`. Defaults to the maximum log likelihood among the `K` clusters using `beta_init`, `tau_init`, `beta_unexch_init`, and `tau_unexch_init`. Note: if `z0` is specified then `eps0` must also be specified.
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
#' ## Fit models stratified by arm (no borrowing)
#' fit.stratified = lapply(
#'   0:1, function(arm) {
#'     tte_reg_bnp(
#'       survival::Surv(y, event) ~ x
#'       , data = curdata[curdata$a == arm, ]
#'       , external_data = histdata[histdata$a == arm, ]
#'       , K = 5, K_unexch = 5
#'       , nburnin = 0, nsamples = 100
#'     )
#'   }
#' )
#' 
#' ## Marginalize strata only (no Bayesian bootstrap)
#' marg.strataonly = tte_marginalize(
#'   fit.strataonly, curdata, stratum_var = 'a'
#'   , hazard_times = times, log_hazard = TRUE
#'   , survival_times = times, log_survival = FALSE
#' )
#' 
#' ## Marginalize models stratified by arm (does Bayesian bootstrap). Note that
#' ## `object` can be a list.
#' marg.stratified = tte_marginalize(
#'   fit.anova_ddp, curdata, stratum_var = 'a'
#'   , hazard_times = times, log_hazard = TRUE
#'   , survival_times = times, log_survival = FALSE
#' )
#' 
#' ## Summarize marginal estimates
#' summary( marg.strataonly )
#' summary( marg.stratified, mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)) )
#' 
#' @export
#' @method summary bnpdb_tte_marginal
summary.bnpdb_tte_marginal <- function(object, ...) {
  Q          = length(object)      ## number of summary quantities
  Qnames     = names(object)       ## names of the summary quantities
  out        = vector('list', Q)
  names(out) = Qnames
  attrib     = attributes(object)
  for ( q in 1:Q ) {
    draws     = object[[q]]
    draws_dim = dim(draws)
    timename  = paste0(Qnames[q], '_times')
    times     = NULL
    if ( timename %in% names(attrib) ) 
      times = attrib[[timename]]
    
    ## Skip summary if there are no draws
    if ( any( draws_dim == 0 ) )
      next
    
    ## Otherwise, perform summary by stratum
    S               = dim(draws)[3]
    Snames          = dimnames(draws)[[3]]
    out[[q]]        = vector('list', S)
    names(out[[q]]) = Snames
    for ( s in 1:S ) {
      tmp           = cbind(summarize_draws(draws[, , s], ...), time = times)
      out[[q]][[s]] <- tmp[, c(1, ncol(tmp), 2:(ncol(tmp) - 1))]
    }
  }
  class(out) = c('bnpdb_tte_marginal_summary', class(out))
  out
}



#' @title Print Method for 'bnpdb_tte_marginal_summary' Objects
#'
#' @description
#' Prints a formatted summary of marginal posterior estimates (e.g., hazard and survival functions)
#' stratified by component and subgroup. Typically called after using \code{summary()} on an object
#' of class \code{bnpdb_tte_marginal}.
#'
#' @param x An object of class \code{bnpdb_tte_marginal_summary}, usually returned by
#'   \code{\link{summary.bnpdb_tte_marginal}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.bnpdb_tte_marginal}}, \code{\link{tte_marginalize}}
#' 
#' @importFrom cli col_blue
#'
#' @keywords internal
#' @export
#' @method print bnpdb_tte_marginal_summary
print.bnpdb_tte_marginal_summary <- function(x, ...) {
  for (quantity_name in names(x)) {
    quantity <- x[[quantity_name]]
    
    # Skip if the quantity is NULL
    if (is.null(quantity)) next
    
    # Loop through strata
    for (stratum_name in names(quantity)) {
      draw_summary <- quantity[[stratum_name]]
      
      # Skip if the element is NULL or has zero rows or columns
      if (is.null(draw_summary) || length(dim(draw_summary)) != 2 || any(dim(draw_summary) == 0)) next
      
      # Header: quantity name + stratum name
      cat(cli::col_blue(paste0(quantity_name, " (", stratum_name, ")")), "\n")
      
      # Print the draw summary table
      print(draw_summary)
      cat("\n")
    }
  }
  invisible(x)
}




#' @title Print Method for 'bnpdb_tte_conditional_summary' Objects
#'
#' @description
#' Prints summary statistics of posterior draws for parameters from a conditional
#' Bayesian nonparametric time-to-event model. Relies on the default method for
#' objects of class \code{posterior::draws_summary}.
#'
#' @param x An object of class \code{bnpdb_tte_conditional_summary}, typically returned by
#'   \code{\link{summary.bnpdb_tte_conditional}}.
#' @param ... Additional arguments passed to the print method.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.bnpdb_tte_conditional}}, \code{\link[posterior]{print.draws_summary}}
#'
#' @keywords internal
#' @export
#' @method print bnpdb_tte_conditional_summary
print.bnpdb_tte_conditional_summary <- function(x, ...) {
  NextMethod("print")
}


#' @title Print Method for 'bnpdb_tte_conditional' Objects
#'
#' @description
#' Prints a summary of a Bayesian nonparametric time-to-event model with conditional borrowing.
#' Internally calls \code{summary()} on the object before printing.
#'
#' @param x An object of class \code{bnpdb_tte_conditional}, as returned by \code{\link{tte_reg_bnp}}.
#' @param ... Additional arguments passed to \code{summary()}.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.bnpdb_tte_conditional}}, \code{\link{tte_reg_bnp}}
#'
#' @keywords internal
#' @export
#' @method print bnpdb_tte_conditional
print.bnpdb_tte_conditional <- function(x, ...) {
  print(summary(x, ...))
}

#' @title Print Method for 'bnpdb_tte_marginal' Objects
#'
#' @description
#' Prints a summary of marginal posterior draws (e.g., survival or hazard functions)
#' from a Bayesian nonparametric time-to-event model. Internally calls \code{summary()} on the object.
#'
#' @param x An object of class \code{bnpdb_tte_marginal}, typically returned by \code{\link{tte_marginalize}}.
#' @param ... Additional arguments passed to \code{summary()}.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.bnpdb_tte_marginal}}, \code{\link{tte_marginalize}}
#'
#' @keywords internal
#' @export
#' @method print bnpdb_tte_marginal
print.bnpdb_tte_marginal <- function(x, ...) {
  print(summary(x, ...))
}


