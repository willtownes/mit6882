sample_hmm_p0 <- 
  function(
    # data & likelihood
    x, log_lik_func, 
    # CRF parameter, prior
    gamma, alpha, kappa, lambda,
    # CRF parameter, data
    z, m, p0, pk, n_jk,
    verbose = FALSE)
    ######################################################
    # sample pi, the global dish distribution
    # input:
    # > lik_func: likelihood function return lik and n_ji
    # > z:  dish assignment for obs
    # > m:  table count per parameter
    # output:
    # > pi: updated z assigment
    ######################################################    
  {
    K <- length(m)
    
    pi_new <- rdirichlet(m + gamma/K)
    
    pi_new
  }