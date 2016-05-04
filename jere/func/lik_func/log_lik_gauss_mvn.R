library(mvtnorm)


log_lik_gauss_mvn <- 
  function(x, lambda){
    # Input
    # x:      data
    # lambda: list of emission parameters (mu, cov)
    
    dmvnorm(x, lambda$mean, lambda$cov, log = TRUE)
  }
