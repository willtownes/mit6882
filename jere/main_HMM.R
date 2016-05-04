require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func/")

set.seed(10)
dat <- 
  HMM_gauss(# 3 cluster, 1 sub-cluster (per cl), 
    gamma = 1, alpha = 10, kappa = 10, 
    K = 10, T = 1000,
    d_gauss = 2, var_gauss = 0.01)

plot(dat$X[, 1])
dat$Z

# run HMM assume true parameters known
# warning: message explode if too many obs on one mode

x <- dat$X
log_lik_func <- log_lik_gauss_mvn
lambda_init <- dat$Theta 

z <- HMM_HDP(
  # data, restaurant number, dish number
  x, K = 10,
  # base measure & parameters
  log_lik_func = log_lik_gauss_mvn, 
  lambda_init = lambda_init, 
  sample_lambda = FALSE, lambda_sampler = NULL, 
  # DP concentration/base-measure
  gamma = 1, alpha = NULL, kappa = NULL, 
  # initialization
  z_init = NULL, m_init = NULL, 
  p0_init = NULL, pk_init = NULL,
  # MCMC parameter
  iter_max = 2e3, iter_update = 10, 
  verbose = FALSE)
