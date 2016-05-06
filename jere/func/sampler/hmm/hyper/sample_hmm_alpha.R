sample_hmm_alpha <- 
  function(
    y, log_lik_func, 
    z, m, x,
    hyper, theta, lambda,
    p0, pk)
  {
    gamma <- hyper$gamma
    alpha <- hyper$alpha
    kappa <- hyper$kappa
    
    #
    alpha_new <- alpha
    
    
    #
    hyper$alpha <- alpha_new
    hyper    
  }