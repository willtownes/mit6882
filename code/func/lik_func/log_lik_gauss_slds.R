library(mvtnorm)


log_lik_gauss_slds <- 
  function(y, x, theta){
    # INPUT
    # y:      observation to be evaluated
    # x:      latent observation to calculate mean
    # lambda: list of emission parameters (A, Sigma)
    dmvnorm(y,
            mean  = theta$A %*% x + theta$B, 
            sigma = theta$Sigma, 
            log = TRUE)
  }
