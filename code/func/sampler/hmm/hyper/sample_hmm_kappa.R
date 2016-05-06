sample_hmm_kappa <- 
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
    kappa_new <- kappa
    
    
    #
    hyper$kappa <- kappa_new
    hyper
  }