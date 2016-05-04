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

plot(dat$X[, 1], col = dat$Z, pch = 19)
z_true <- dat$Z

# run HMM assume true parameters known
# warning: message explode if too many obs on one mode

y <- dat$X
log_lik_func <- log_lik_gauss_mvn
theta_init <- dat$Theta 
lambda_init <- 0

theta_sampler <- sample_theta
lambda_sampler <- sample_lambda
x_sampler <- sample_x

z_obs <- 
  HMM_HDP(
    # data, restaurant number, dish number
    y, K = 10,
    # emission measure & parameters
    log_lik_func = log_lik_func, 
    sample_emiss = FALSE, 
    theta_init = theta_init, 
    theta_sampler = theta_sampler, 
    lambda_init = lambda_init, 
    lambda_sampler = lambda_sampler,
    x_sampler = x_sampler,
    # DP concentration/base-measure
    sample_hyper = FALSE,
    gamma = 1, alpha = NULL, kappa = NULL, 
    # initialization
    x_init = NULL, z_init = NULL, m_init = NULL, 
    p0_init = NULL, pk_init = NULL,
    # MCMC parameter
    iter_max = 10, iter_update = 1, verbose = TRUE)

table(as.factor(z_obs[9, ]), as.factor(z_true))
