require(magrittr)
require(dplyr)
source("./func/util/source_Dir.R")
sourceDir("./func/")

load("will_dat.RData")
pars[[1]]$B <- pars[[1]]$B/100
pars[[1]]$Sigma <- pars[[1]]$Sigma/1000
pars[[2]]$A <- pars[[2]]$A
pars[[2]]$B <- pars[[2]]$B/500
pars[[2]]$Sigma <- pars[[2]]$Sigma/3000


set.seed(8000)
z_true <- 
  HMM_z(gamma = 1, alpha = 1, kappa = 10, 
        K = 2, T = 100)

z_true <- c(rep(1, 30), rep(2, 30), rep(1, 30), rep(2, 30))
dat <- HMM_SLDS(z_true, pars)

plot(dat$Y, col = dat$Z, pch = 19)
z_true <- dat$Z

# generate y


# run HMM assume true parameters known
# warning: message explode if too many obs on one mode

y <- dat$Y
K <- 2
x_init <- dat$X
x_sampler <- sample_x

log_lik_func <- log_lik_gauss_slds

theta_init <- dat$Theta 
lambda_init <- 0

theta_sampler <- sample_theta
lambda_sampler <- sample_lambda

z_obs <- 
  HMM_HDP(
    # data, restaurant number, dish number
    y, K = 2,
    # emission measure & parameters
    log_lik_func = log_lik_func, 
    sample_emiss = FALSE, 
    theta_init = theta_init, theta_sampler = theta_sampler, 
    lambda_init = lambda_init, lambda_sampler = lambda_sampler,
    x_init = x_init, x_sampler = x_sampler,
    # DP concentration/base-measure
    sample_hyper = FALSE,
    gamma = 1, alpha = NULL, kappa = NULL, 
    # initialization
    z_init = NULL, m_init = NULL, 
    p0_init = NULL, pk_init = NULL,
    # MCMC parameter
    iter_max = 10, iter_update = 1, verbose = TRUE)

table(as.factor(z_obs[nrow(z_obs), ]), 
      as.factor(z_true))
