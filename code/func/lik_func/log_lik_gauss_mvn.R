library(mvtnorm)


log_lik_gauss_mvn <- 
  function(y, theta){
    # Input
    # x:      data
    # lambda: list of emission parameters (mu, cov)
    
    dmvnorm(y, theta$mean, theta$cov, log = TRUE)
  }
