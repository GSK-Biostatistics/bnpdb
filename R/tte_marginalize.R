
#' Function to obtain RHS formula in case user omits outcome in data
#' @keywords internal
#' @noRd
get_rhs_formula <- function(formula) {
  term_labels <- attr(terms(formula), "term.labels")
  if (length(term_labels) == 0) {
    return(~ 1)
  } else {
    return(reformulate(term_labels))
  }
}

#' Formula to make design matrices based on fixing a certain strata to be one stratum
#' @keywords internal
#' @noRd
make_stratified_design_matrices <- function(data, formula, stratum_var) {

  ## Check if stratum_var refers to a variable that is coded as 0:1 or
  ## a factor. The following are acceptable:
  ##   y ~ a + ...;         a is numeric or integer with range 0:1
  ##   y ~ a + ...          a is coded as a factor variable
  ##   y ~ factor(a) + ...  a is anything
  formula = get_rhs_formula(formula)
  mf      = model.frame(formula, data = data)
  terms   = terms(mf)
  dc      = attr(terms, 'dataClasses')
  xvars   = all.vars(formula)
  strata  = data[[stratum_var]]
  
  ## If stratum variable is not in the formula, return the regular design matrix
  ## as a list
  if ( !( stratum_var %in% xvars ) ) {
    return( lapply() )
  } else {  
    ## Otherwise, must create list of stratum-specific design matrices where the
    ## stratum membership is fixed
    
  }
  
  ## Otherwise, must create 
  stratum_var.class = dc[xvars == stratum_var]  ## strata data type as determined by data and formula

  if ( stratum_var.class %in% c('numeric', 'integer') ) {
    ## Throw an error if the stratum doesn't make sense.
    stopifnot(
        "`stratum_var` refers to a numeric variable outside the range 0:1."
      = all(strata %in% c(0, 1))
    )

    ## Recode to be a factor so that the rest of the code generalizes the
    ## result
    data[[stratum_var]] = factor(data[[stratum_var]], levels = c(0, 1))
  }


  ## It must now be the case that one of the following is true:
  ##     y ~ a + ...          a is coded as a factor variable
  ##     y ~ factor(a) + ...  a is anything (numeric, etc.)
  ## In either case, we are working with a factor variable
  ## so we convert the data to reflect this.
  data[[stratum_var]] = factor(data[[stratum_var]])

  ## Get number of strata and repeat data that many times
  levs    = levels(data[[stratum_var]])
  S       = length(levs)
  n       = nrow(data)
  repdata = lapply(levs, function(x) {
    sdata = data
    sdata[[stratum_var]] = x
    return(sdata)
  })
  repdata = do.call('rbind', repdata)

  ## Now, obtain design matrix based on the stacked repdata
  Xlong = model.matrix(formula, repdata)

  ## Convert Xlong to a list of design matrices: one per stratum
  Xlist = lapply(
    1:S, function(s) {
      end   = s * n
      start = end - n + 1
      Xlong[start:end, ]
    }
  )
  Xlist
}


#' Marginalize DPMM output for time-to-event data
#' 
#' Report posterior hazard, survival, and/or RMST functions based on Dirichlet
#' process mixture model output. Uses the Bayesain bootstrap when there are
#' covariates
#' 
#' @param samples output from calling \code{\link{tte_reg_bnp}}
#' @param data a \code{\link[base]{data.frame}} consisting of covariates over which marginalization should be performed
#' @param stratum_var a string giving the name of the stratification variable (e.g., treatment arm) present in the `formula` used when calling \code{\link{tte_reg_bnp}} as well as the `data`
#' @param hazard_times a \code{\link[base]{numeric}} `vector` of nonnegative values giving the time points to evaluate the hazard function; if `NULL`, will not evaluate hazard function
#' @param log_hazard `logical`: if `TRUE`, evaluates the hazard function on the log scale; ignored if `hazard_times` is NULL
#' @param survival_times a \code{\link[base]{numeric}} `vector` of nonnegative values giving the time points to evaluate the survival function; if `NULL`, will not evaluate survival function
#' @param log_survival `logical`: if `TRUE`, evaluates the survival function on the log scale; ignored if `survival_times` is NULL
#' @param rmst_times a \code{\link[base]{numeric}} `vector` of nonnegative values giving the time points to evaluate the restricted mean survival time (RMST); if `NULL`, will not evaluate the RMST
#' @param log_rmst `logical`: if `TRUE`, evaluates the RMST on the log scale; ignored if `rmst_times` is NULL
#' 
#' @examples
#' 
#'  ##  Set seed for reproducibility
#'  set.seed(741)
#'  
#'  ##  Generate current and historical data sets
#'  n     <- 100
#'  n0    <- 100
#'  N     <- n + n0
#'  pexch <- 0.80
#'  beta  <- cbind(c(1, 1, 0), c(-1, -1, 0))
#'  sigma <- c(0.5, 0.75)
#'  a     <- rbinom(N, 1, 0.5)      ##  tretment indicator
#'  x     <- rbinom(N, 1, 0.5)      ##  binary covariate
#'  eps0  <- rbinom(n0, 1, pexch)   ##  exchangeability indicator
#'  eps   <- c(rep(1, n), 2 - eps0) ##  1 = exch; 2 = unexch
#'  X     <- cbind(1, a, x)
#'  Mu    <- unlist( 
#'    lapply(1:N, function(i) { as.numeric(X[i, ] %*% beta[, eps[i]] ) } ) 
#'  )
#'  logt     <- rnorm(N, Mu, sigma[eps])
#'  logc     <- rnorm(N, mean(Mu) + 0.50, sd = 0.25)
#'  event    <- logt <= logc
#'  logy     <- ifelse(event, logt, logc)
#'  dat      <- data.frame(y = exp(logy), event = event, a = a, x = x)
#'  curdata  <- dat[1:n, ]
#'  histdata <- dat[-(1:n), ]
#'  
#'  ##  Times for marginalization
#'  tau   = max(curdata$y[curdata$event == 1])
#'  times = seq(0.01, tau, length.out = 5)
#'  
#'  ##  Fit intercept only (no borrowing)
#'  fit.intonly = tte_reg_bnp(
#'    survival::Surv(y, event) ~ 1, data = curdata, K = 5
#'    , nburnin = 0, nsamples = 100
#'  )
#'  
#'  ##  Fit treatment only (with borrowing)
#'  fit.strataonly = tte_reg_bnp(
#'    survival::Surv(y, event) ~ a, data = curdata, K = 5
#'    , external_data = histdata, K_unexch = 5
#'    , nburnin = 0, nsamples = 100
#'  )
#'  
#'  
#'  ##  Fit treatment and binary covariate (with borrowing)
#'  fit.anova_ddp = tte_reg_bnp(
#'    survival::Surv(y, event) ~ a + x, data = curdata, K = 5
#'    , external_data = histdata, K_unexch = 5
#'    , nburnin = 0, nsamples = 100
#'  )
#'  
#'  ##  Fit models stratified by arm (no borrowing)
#'  fit.stratified = lapply(
#'    0:1, function(arm) {
#'      tte_reg_bnp(
#'        survival::Surv(y, event) ~ x
#'        , data = curdata[curdata$a == arm, ]
#'        , external_data = histdata[histdata$a == arm, ]
#'        , K = 5, K_unexch = 5
#'        , nburnin = 0, nsamples = 100
#'      )
#'    }
#'  )
#'  
#'  ##  Marginalize strata only (no Bayesian bootstrap)
#'  marg.strataonly = tte_marginalize(
#'    fit.strataonly, curdata, stratum_var = 'a'
#'    , hazard_times = times, log_hazard = TRUE
#'    , survival_times = times, log_survival = FALSE
#'  )
#'  
#'  ##  Marginalize w/ covariate model(uses Bayesian bootstrap)
#'  marg.anova_ddp = tte_marginalize(
#'    fit.anova_ddp, curdata, stratum_var = 'a'
#'    , hazard_times = times, log_hazard = TRUE
#'    , survival_times = times, log_survival = FALSE
#'  )
#'  
#'  ## Marginalize with stratified covariate model (uses Bayesian bootstrap)
#'  marg.stratified = tte_marginalize(
#'    fit.stratified, curdata, stratum_var = 'a'
#'    , hazard_times = times, log_hazard = TRUE
#'    , survival_times = times, log_survival = FALSE
#'  )
#'  
#'  ##  Summarize marginal estimates
#'  summary( marg.strataonly )
#'  summary( marg.anova_ddp, mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)) )
#'  summary( marg.stratified, mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.975)) )
#'
#' @return A named list of class \code{tte_marginalize} consisting of marginalized posterior samples.
#' Let \eqn{M}, \eqn{T}, and \eqn{S} denote the number of MCMC samples, time points, and strata, respectively.
#' The returned object contains:
#' \itemize{
#'   \item \code{survival}: an \eqn{M \times T \times S} array giving posterior samples for the survival function for each stratum
#'   \item \code{hazard}: an \eqn{M \times T \times S} array giving posterior samples for the hazard function for each stratum
#'   \item \code{rmst}: an \eqn{M \times T \times S} array giving posterior samples for the restricted mean survival time (RMST) for each stratum
#' }
#' 
#' @export
tte_marginalize <- function(
    samples, data, stratum_var = NULL
    , hazard_times = NULL, log_hazard = TRUE
    , survival_times = NULL, log_survival = TRUE
    , rmst_times = NULL, log_rmst = TRUE
  ) {
  
  ## Sort data by stratum_var
  if ( !is.null(stratum_var) )
    data = data[order(data[[stratum_var]]), ]
  
  ## Check inputs
  if ( is.null(hazard_times) & is.null(survival_times) & is.null(rmst_times) )
    stop("At least one of `hazard_times`, `survival_times`, or `rmst_times` must be specified.")
  
  stopifnot("`samples` should be an output from `tte_reg_bnp` or a list of outputs from `tte_reg_bnp`" = is.list(samples))
  
  if ( !is.null(hazard_times) ) {
    stopifnot("`hazard_times` should be a vector of nonnegative value if specified" = { length(hazard_times) > 1; 'numeric' %in% class(hazard_times); complete.cases(hazard_times) })
  } else hazard_times = numeric(0)
  
  if ( !is.null(survival_times) ) {
    stopifnot("`survival_times` should be a vector of nonnegative value if specified" = { length(survival_times) > 1; 'numeric' %in% class(survival_times); complete.cases(survival_times) })
  } else { survival_times = numeric(0) }
  
  if ( !is.null(rmst_times) ) {
    stopifnot("`rmst_times` should be a vector of nonnegative value if specified" = { !is.null(rmst_times); length(rmst_times) > 1; 'numeric' %in% class(rmst_times); complete.cases(rmst_times) })
  } else rmst_times = numeric(0)
  
  ## Check stratum_var
  if( !is.null(stratum_var) ) {
    stopifnot("If specified, `stratum_var` must be a character of length 1 giving a variable in `data` that contains the name of the stratification factor (e.g, treatment arm)" = { is.character(stratum_var); length(stratum_var) == 1; stratum_var %in% names(data) })
  } else {
    stopifnot( "If `stratum_var` is NULL, `samples` must refer to only a single posterior" = ( 'bnpdb_tte_conditional' %in% class(samples) ) )
  }

  
  ## Check if there is only one set of posterior samples
  posterior_single = FALSE
  if ( all( c('beta', 'sigma', 'alpha', 'w') %in% names(samples) ) ) {
    posterior_single = TRUE
    samples = list(samples)
  }
  if ( length(samples) == 1 & 'bnpdb' %in% class(samples[[1]]) ) {
    posterior_single = TRUE
  }
  
  ## Obtain all formulas
  rhs.formula.list = lapply(samples, function(x) get_rhs_formula( attr(x, 'formula')) )
  
  
  ## Check consistency of the formulas--all must either have stratum in the
  ## formula or not a single one can.
  has_stratum <- vapply(
    rhs.formula.list,
    function(f) {
      vars <- all.vars(f)
      if (length(vars) == 0) {
        FALSE
      } else {
        stratum_var %in% vars
      }
    },
    logical(1)
  )
  
  if (any(has_stratum) && !all(has_stratum)) {
    stop("Inconsistent use of `stratum_var` across models: it must be either included in all or excluded from all.")
  }
  
  ## Check if Bayes bootstrap is necessary 
  ## (it is if the formula includes anything other than the stratum variable / intercept)
  bayesBoot = sapply(rhs.formula.list, function(f) { 
    rhs_vars = all.vars(f)
    !(length(rhs_vars) == 0 || (length(rhs_vars) == 1 && rhs_vars == stratum_var))
  } )
  bayesBoot = any(bayesBoot)   ## convert to 1 dimension
  
  
  
  ## If Bayesian bootstrap is not necessary, call the C++ function that uses
  ## a single design matrix and does not sample the bootstrap weights.
  
  if ( !bayesBoot ) {  ## Stratum-only or intercept-only model
    
    ## Get unique design matrix--sort by ascending stratum variable
    X  = unique( model.matrix(rhs.formula.list[[1]], data = data) )
    
    ## Call function based on X as a matrix
    out = marginalize_strataonly(
      X, samples[[1]]
      , log(survival_times), log_survival
      , log(hazard_times), log_hazard
      , log(rmst_times), log_rmst
    )
    
    
  } else {  ## Need to integrate over covariate distribution
    
    ## Obtain list of design matrices: one for each estimated model / stratum.
    if ( all(has_stratum) ) { 
      ## if the stratum variable is in the formula:
      Xlist   = make_stratified_design_matrices(data, rhs.formula.list[[1]], stratum_var)
      samples = lapply(1:length(unique(data[[stratum_var]])), function(i) samples[[1]])
    } else {
      Xlist = lapply(rhs.formula.list, function(f) model.matrix(f, data))
    }
    
    ## Call g-computation C++ function to perform marginalization
    out = marginalize_tte_gcomp(
      Xlist, samples
      , log(survival_times), log_survival
      , log(hazard_times), log_hazard
      , log(rmst_times), log_rmst
    )
  }
  
  ## Obtain stratum names for output
  if ( !is.null(stratum_var) ) {
    snames = paste0(stratum_var, ' = ', unique(data[[stratum_var]]) )
  } else {
    snames = ""
  }
  
  ## Dimension names
  if ( length(hazard_times) > 0 ) {
    prefix = ifelse(log_hazard, 'log_haz[', 'haz[')
    colnames(out[['hazard']])         = paste0(prefix, 1:length(hazard_times), ']')
    dimnames(out[['hazard']])[3][[1]] = snames
  } else { out[['hazard']] = array(dim = c(0, 0, 0)) }
  
  if ( length(survival_times) > 0 ) {
    prefix = ifelse(log_survival, 'log_surv[', 'surv[')
    colnames(out[['survival']])         = paste0(prefix, 1:length(survival_times), ']')
    dimnames(out[['survival']])[3][[1]] = snames
  } else { out[['survival']] = array(dim = c(0, 0, 0)) }
  
  if ( length(rmst_times) > 0 ) {
    prefix = ifelse(log_rmst, 'log_rmst[', 'rmst[')
    colnames(out[['rmst']])         = paste0(prefix, 1:length(rmst_times), ']')
    dimnames(out[['rmst']])[3][[1]] = snames
  } else { out[['rmst']] = array(dim = c(0, 0, 0)) }
  
  ## Change class and enrich with attributes for some possible downstream
  ## analyses
  class(out) = c('bnpdb', class(out))
  attr(out, "datatype")       = "tte"
  attr(out, "type")           = "marginal parameters"
  attr(out, "stratum_var")    = stratum_var
  attr(out, "hazard_times")   = hazard_times
  attr(out, "log_hazard")     = log_hazard
  attr(out, "survival_times") = survival_times
  attr(out, "log_survival")   = log_survival
  attr(out, "rmst_times")     = rmst_times
  attr(out, "log_rmst")       = log_rmst
  class(out) = c("bnpdb_tte_marginal", class(out))

  return(out)
}
