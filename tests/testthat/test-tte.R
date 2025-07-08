##  Generate current and historical data sets
set.seed(741)
n     = 100
n0    = 100
N     = n + n0
pexch = 0.80
beta  = cbind(c(1, 1, 0), c(-1, -1, 0))
sigma = c(0.5, 0.75)
a     = rbinom(N, 1, 0.5)      ##  tretment indicator
x     = rbinom(N, 1, 0.5)      ##  binary covariate
eps0  = rbinom(n0, 1, pexch)   ##  exchangeability indicator
eps   = c(rep(1, n), 2 - eps0) ##  1 = exch; 2 = unexch
X     = cbind(1, a, x)
Mu    = unlist( 
  lapply(1:N, function(i) { as.numeric(X[i, ] %*% beta[, eps[i]] ) } ) 
)
logt     = rnorm(N, Mu, sigma[eps])
logc     = rnorm(N, mean(Mu) + 0.50, sd = 0.25)
event    = logt <= logc
logy     = ifelse(event, logt, logc)
dat      = data.frame(y = exp(logy), event = event, a = a, x = x)
curdata  = dat[1:n, ]
histdata = dat[-(1:n), ]
nsamples = 200

##  Times for marginalization
tau   = max(curdata$y[curdata$event == 1])
times = seq(0.01, tau, length.out = 5)


test_that("DPMM returns posterior samples", {
  ##  Fit intercept only (no borrowing)
  fit.intonly = tte_reg_bnp(
    Surv(y, event) ~ 1, data = curdata, K = 5
    , nburnin = 0, nsamples = nsamples
  )
  
  ##  Fit treatment only (with borrowing)
  fit.strataonly = tte_reg_bnp(
    Surv(y, event) ~ a, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  
  
  ##  Fit treatment and binary covariate (with borrowing)
  fit.anova_ddp = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  
  ##  Fit models stratified by arm (no borrowing)
  fit.stratified = lapply(
    0:1, function(arm) {
      tte_reg_bnp(
        Surv(y, event) ~ x
        , data = curdata[curdata$a == arm, ]
        , external_data = histdata[histdata$a == arm, ]
        , K = 5, K_unexch = 5
        , nburnin = 0, nsamples = nsamples
      )
    }
  )
  expect_type(fit.intonly, "list")
  exchpar_names   = c('beta', 'sigma', 'alpha', 'w')
  unexchpar_names = c('eps0', 'pexch')
  expect_true( all( exchpar_names %in% names(fit.intonly) ) )
  expect_true( all( c(exchpar_names, unexchpar_names) %in% names(fit.strataonly) ) )
  expect_true( all( c(exchpar_names, unexchpar_names) %in% names(fit.anova_ddp) ) )
  expect_true( ncol(fit.strataonly$eps0) == n0 )
})

test_that("External data labels are integers", {
  fit = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  expect_true(is.numeric(fit$eps0))
  expect_true(all(fit$eps0 == floor(fit$eps0)))
})

test_that("Reproducibility with fixed seed", {
  set.seed(42)
  fit1 = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  
  set.seed(42)
  fit2 = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  
  expect_equal(fit1$clusters, fit2$clusters)
})

test_that("DPMM handles empty input gracefully", {
  expect_error(tte_reg_bnp(NULL))
})

test_that("Conditional summary functions working", {
  fit = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  fit.summary = summary(fit, mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)))
  expect_true( 'bnpdb_tte_conditional_summary' %in% class(fit.summary) )
  expect_true( 'draws_summary' %in% class(fit.summary) )
})


test_that("Marginalization without Bayesian bootstrap", {
  fit = tte_reg_bnp(
    Surv(y, event) ~ a, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  marg = tte_marginalize(
    fit, data = curdata, stratum_var = 'a'
    , hazard_times = times, log_hazard = TRUE
    , survival_times = times, log_survival = FALSE
  )
  summ = summary(marg)
  
  expect_true("bnpdb_tte_marginal" %in% class(marg))
  expect_true( all( c("survival", "rmst", "hazard") %in% names(marg) ) )
  expect_true( all( dim(marg$hazard) == c(nsamples, 5, 2) ) )
  expect_true( class(summ$survival) == 'list' )
  expect_true( all( summ$hazard$`a = 0`$variable == paste0('log_haz[', 1:5, ']') ) )
  expect_true( all( summ$survival$`a = 0`$variable == paste0('surv[', 1:5, ']') ) )
})


test_that("Marginalization with Bayesian bootstrap", {
  fit = tte_reg_bnp(
    Surv(y, event) ~ a + x, data = curdata, K = 5
    , external_data = histdata, K_unexch = 5
    , nburnin = 0, nsamples = nsamples
  )
  marg = tte_marginalize(
    fit, data = curdata, stratum_var = 'a'
    , hazard_times = times, log_hazard = TRUE
    , survival_times = times, log_survival = FALSE
  )
  summ = summary(marg)
  
  expect_true("bnpdb_tte_marginal" %in% class(marg))
  expect_true( all( c("survival", "rmst", "hazard") %in% names(marg) ) )
  expect_true( all( dim(marg$hazard) == c(nsamples, 5, 2) ) )
  expect_true( class(summ$survival) == 'list' )
  expect_true( all( summ$hazard$`a = 0`$variable == paste0('log_haz[', 1:5, ']') ) )
  expect_true( all( summ$survival$`a = 0`$variable == paste0('surv[', 1:5, ']') ) )
})
