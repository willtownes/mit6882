sample_hdp_pi <- 
  function(
    # data & likelihood
    x, lik_func, 
    # CRF parameter, prior
    gamma, alpha, J, 
    # CRF parameter, data
    z, m, pi)
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
    pi_new <- rdirichlet(c(m, gamma))
    pi_new
  }