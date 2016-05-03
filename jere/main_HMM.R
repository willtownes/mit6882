require(magrittr)
require(dplyr)
load("will_data.RData")

source("./func/util/source_Dir.R")
sourceDir("./func/")
lik_func <- lik_gauss_mean

set.seed(10)
dat <- 
  HMM_gauss(# 3 cluster, 1 sub-cluster (per cl), 
    gamma = 1, alpha = 10, kappa = 10, 
    K = 10, T = 100,
    var_gauss = 0.01)

plot(dat$X, type = "l")


HDP(
  # data, restaurant number, dish number
  x, nJ = 1, nK_init = 1,
  # DP concentration/base-measure
  gamma = 0.1, alpha = 0.5, lik_func,
  # initialization
  z_init = NULL, m_init = NULL, pi_init = NULL,
  # MCMC parameter
  iter_max = 1e3,
  # utility
  J_init = c("random", "kmeans")[1]
)
