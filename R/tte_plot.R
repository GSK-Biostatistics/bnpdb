

#' Plot time-to-event posterior
#'
#' Plot the posterior for time-to-event data
#'
#' @param marginal_samples output from calling \code{\link{tte_marginalize}}
#' @param point_estimate quantity to use as the point estimate when plotting
#' @param ci_prob credible interval bounds to plot; if `NULL`, only plots the point estimate
#' @param quantity a `character` giving which posterior quantities to plot
#' 
#' @import ggplot2
#' @importFrom dplyr mutate arrange group_by
#' 
#' @return a list of posterior function plots
#' @examples
#' ##  Generate current and historical data sets
#' set.seed(741)
#' n     <- 100
#' n0    <- 100
#' N     <- n + n0
#' pexch <- 0.80
#' beta  <- cbind(c(1, 1, 0), c(-1, -1, 0))
#' sigma <- c(0.5, 0.75)
#' a     <- rbinom(N, 1, 0.5)      ##  tretment indicator
#' x     <- rbinom(N, 1, 0.5)      ##  binary covariate
#' eps0  <- rbinom(n0, 1, pexch)   ##  exchangeability indicator
#' eps   <- c(rep(1, n), 2 - eps0) ##  1 = exch; 2 = unexch
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
#' ##  Times for marginalization
#' tau   = max(curdata$y[curdata$event == 1])
#' times = seq(0.01, tau, length.out = 5)
#' 
#' ##  Fit treatment only (with borrowing)
#' fit.strataonly = tte_reg_bnp(
#'   survival::Surv(y, event) ~ a, data = curdata, K = 5
#'   , external_data = histdata, K_unexch = 5
#'   , nburnin = 0, nsamples = 100
#' )
#' 
#' ##  Marginalize strata only (no Bayesian bootstrap)
#' marg.strataonly = tte_marginalize(
#'   fit.strataonly, curdata, stratum_var = 'a'
#'   , hazard_times = times, log_hazard = TRUE
#'   , survival_times = times, log_survival = FALSE
#' )
#' 
#' #' plot
#' suppressWarnings({
#'   plotlist = tte_plot(marg.strataonly)
#' })
#' 
#' @export
tte_plot = function(
  marginal_samples, point_estimate = 'mean', ci_prob = 0.95
  , quantity = c('survival', 'hazard', 'rmst')
) {
  stopifnot("`marginal_samples` must be an output from calling `tte_marginalize`" = "bnpdb_tte_marginal" %in% class(marginal_samples))
  if (!is.null(ci_prob) ) {
    stopifnot("`ci_prob` must be NULL or a number larger than 0 and less than 1" = { ci_prob > 0; ci_prob < 1 } )
  }
  attrib    = attributes(marginal_samples)
  stopifnot( "`quantity` must be at least one of `survival`, `hazard` or `rmst`" = all( quantity %in% c('survival', 'hazard', 'rmst')) )

  ## Call summary function
  if ( !is.null(ci_prob) ) {
    summ  = summary(
      marginal_samples, point_estimate
      , lower = ~posterior::quantile2(.x, probs = (1 - ci_prob)/2, names = FALSE )
      , upper = ~posterior::quantile2(.x, probs = 1 - (1 - ci_prob)/2, names = FALSE )
    )
  } else {
    summ  = summary(marginal_samples, point_estimate)
  }
  
  has_samples = sapply(summ, function(x) length(x[[1]]) > 0 )
  to_plot     = quantity[ quantity %in% names(which(has_samples)) ]

  ## Proceed to plot those for which we have samples
  plotlist = list()
  for ( q in to_plot ) {
    quant    = summ[[q]]
    
    ## Check if the summary was on the log scale
    log_ind  = attrib[[paste0('log_', q)]]
    lab      = ifelse(log_ind, paste(q, '(log)'), q )
    
    snames   = names(quant)
    for ( i in 1:length(snames) )
      quant[[i]]$stratum = snames[i]
    plotdata = do.call('rbind', quant)
    names(plotdata)[which( names(plotdata) == point_estimate )] = 'point_estimate'
    
    
    if ( is.null(ci_prob) ) {
      plotlist[[q]] = ggplot2::ggplot(
        data = plotdata
        , mapping = ggplot2::aes(
          x = .data$time, y = .data$point_estimate
          , color = .data$stratum, fill = .data$stratum
        )
      ) +
        ggplot2::geom_point() + 
        ggplot2::geom_smooth(se = FALSE) + 
        ggplot2::ylab(lab) +
        ggplot2::xlab("time") + 
        ggplot2::ggtitle(lab) + 
        theme(legend.title = element_blank())
    } else {
      plotdata = dplyr::mutate(
        dplyr::arrange(
          dplyr::group_by(plotdata, plotdata$stratum),
          time
        ),
        smooth_lower = stats::loess(lower ~ time)$fitted,
        smooth_upper = stats::loess(upper ~ time)$fitted
      )
      
      plotlist[[q]] = ggplot2::ggplot(
        data = plotdata
        , mapping = ggplot2::aes(
          x = .data$time, y = .data$point_estimate
          , color = .data$stratum, fill = .data$stratum
        )
      ) +
        ggplot2::geom_ribbon(
          mapping = ggplot2::aes(ymin = .data$smooth_lower, ymax = .data$smooth_upper),
          alpha = 0.2,
          color = NA
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(se = FALSE) +
        ggplot2::ylab(lab) +
        ggplot2::xlab("time") + 
        ggplot2::ggtitle(lab) + 
        theme(legend.title = element_blank())
    }
  }
  plotlist
}


