sample_hmm_p0 <- 
  function(
    y, log_lik_func, 
    z, m, x,
    p0, pk,
    theta, lambda,
    hyper
  )
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
    gamma <- hyper$gamma
    K <- length(m)
    
    #
    pi_new <- rdirichlet(m + gamma/K)
    
    pi_new
  }