require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func/")

load("will_dat.RData")
tune_will_data <- FALSE
use_will_data <- FALSE
if (tune_will_data){
  theta[[1]]$B <- theta[[1]]$B/50
  theta[[1]]$Sigma <- theta[[1]]$Sigma/1000
  theta[[2]]$A <- theta[[2]]$A
  theta[[2]]$B <- theta[[2]]$B/500
  theta[[2]]$Sigma <- theta[[2]]$Sigma/3000
}

set.seed(100)
# z_true <- 
#   HMM_z(gamma = 1, alpha = 10, kappa = 10, 
#         K = 2, T = 100)
z_true <- c(rep(1, 50), rep(2, 50), rep(1, 50), rep(2, 50), 
            rep(1, 50), rep(2, 50), rep(1, 50), rep(2, 50))
dat <- HMM_SLDS(z_true, theta, x_0 = x_0)
plot(dat$Y, col = dat$Z, pch = 19, cex = 0.5)

plot(dat$Y, pch = 19, cex = 0.5, xlab = "x", ylab = "y")


if (use_will_data){
  z_true <- c(rep(1, 100), rep(2, 100))
  dat <- list(X = as.matrix(rbind(x_0, x_true)), 
              Y = as.matrix(x_true[, 1:2]),
              Theta = theta, Z = z_true)
}

# run HMM assume true parameters known
# warning: message explode if too many obs on one mode
# warning: make sure x include x_0!
y <- dat$Y
K <- 2
x_init <- dat$X
x_sampler<- sample_x

log_lik_func <- log_lik_gauss_slds

theta_init <- dat$Theta 
lambda_init <- 
  list(R = diag(ncol(y)), 
       C = cbind(diag(ncol(y)), 
                 matrix(0, nrow = ncol(y), 
                        ncol = ncol(x_init) - ncol(y))) )

theta_sampler <- sample_theta
lambda_sampler <- sample_lambda

out_list <- 
  HMM_HDP(
    # data, restaurant number, dish number
    y, K = 2,
    # emission measure & parameters
    log_lik_func = log_lik_func, 
    sample_emiss = TRUE, 
    theta_init = theta_init, theta_sampler = theta_sampler, 
    lambda_init = lambda_init, lambda_sampler = lambda_sampler,
    x_init = x_init, x_sampler = x_sampler,
    # DP concentration/base-measure
    sample_hyper = FALSE,
    gamma = 1, alpha = NULL, kappa = 10, 
    # initialization
    z_init = NULL, m_init = NULL, 
    p0_init = NULL, pk_init = NULL,
    # MCMC parameter
    iter_max = 50, iter_update = 1, verbose = FALSE)

z_obs <- out_list$z
table(as.factor(z_obs), as.factor(z_true))

plot(dat$Y[, 1], col = z_obs, pch = 19, cex = 0.5, 
     xlab = "Time", ylab = "x")
plot(dat$Y[, 2], col = z_obs, pch = 19, cex = 0.5,
     xlab = "Time", ylab = "y")


#
z_true <- rep(1, 200)

dat_list <- 
  cbind(runif(100, 0, 10), runif(100, 0, 10)) %>% 
  t %>% as.data.frame %>%
  lapply(function(x_0) HMM_SLDS(z_true, theta, x_0)$X[, 1:2])

plot(0,0, xlim = c(-1e4, 1e4), ylim = c(-1e4, 1e4), 
     xlab = "x", ylab = "y", main = "Example of 'evade' motion")
lapply(dat_list, lines, col = rgb(0,0,0,0.2), lwd = 2)
