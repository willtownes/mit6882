sample_hmm_pk <- 
  function(
    y, log_lik_func, 
    z, m, x,
    p0, pk,
    theta, lambda,
    hyper
  )
    ######################################################
    # sample pk, the state-transition distribution
    ######################################################    
  {
    #
    K <- length(m)
    gamma <- hyper$gamma
    alpha <- hyper$alpha
    kappa <- hyper$kappa
    n_jk <- hyper$n_jk
    
    #
    pk_prob_mat <- alpha * pk + n_jk + kappa * diag(K)
    
    pk_new <- apply(pk_prob_mat, 1, rdirichlet)
    
    t(pk_new)
  }