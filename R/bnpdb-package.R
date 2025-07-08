#' @useDynLib bnpdb, .registration = TRUE
#' @import Rcpp
#' @import survival
#' @import distributional
#' @importFrom cli col_blue
#' @importFrom mvtnorm rmvnorm
#' @importFrom posterior quantile2 summarize_draws
NULL

#' The `bnpdb` Package
#'
#' `bnpdb` provides tools for Bayesian nonparametric modeling of time-to-event data using
#' Dirichlet process mixture models (DPMMs), with support for dynamic borrowing of information 
#' from external data sources. These models allow flexible modeling of complex survival data 
#' with latent subpopulations and nonparametric prior distributions on the space of densities.
#'
#' Features of the package include:
#' \itemize{
#'   \item Nonparametric regression using DPMMs.
#'   \item Posterior summaries for conditional and marginal survival, hazard, and RMST functions.
#'   \item Integration of internal and external data sources with a latent exchangeability framework.
#'   \item Tools for visualizing and comparing model-based inferences.
#' }
#'
#' @references
#' Blackwell, D. and MacQueen, J. (1973). Ferguson distributions via Polya urn schemes. *Annals of Statistics*, **1**, 353–355.  
#' Ferguson, T. S. (1974). Prior distributions on the spaces of probability measures. *Annals of Statistics*, **2**, 615–629.  
#' Sethuraman, J. (1994). A constructive definition of Dirichlet priors. *Statistica Sinica*, **2**, 639–650.  
#' Lo, A. Y. (1984). On a class of Bayesian nonparametric estimates I: Density estimates. *Annals of Statistics*, **12**, 351–357.  
#' Escobar, M. D. (1994). Estimating normal means with a Dirichlet process prior. *Journal of the American Statistical Association*, **89**, 268–277.  
#' Escobar, M. D. and West, M. (1995). Bayesian density estimation and inference using mixtures. *Journal of the American Statistical Association*, **90**, 577–588.  
#' Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models. *Journal of Computational and Graphical Statistics*, **9**, 249–265.  
#' Ishwaran, H. and James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. *Journal of the American Statistical Association*, **96**, 161–173.  
#' Ishwaran, H. and James, L. F. (2002). Approximate Dirichlet process computing in finite normal mixtures: smoothing and prior information. *Journal of Computational and Graphical Statistics*, **11**, 508–532.  
#' Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. A., and Ibrahim, J. G. (2024). LEAP: the latent exchangeability prior for borrowing information from historical data. *Biometrics*, **80**(3), ujae083.  
#'
#' @docType package
#' @name bnpdb
NULL
