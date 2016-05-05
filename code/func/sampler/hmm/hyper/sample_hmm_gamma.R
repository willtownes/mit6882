sample_hmm_gamma <- 
  function(
    y, log_lik_func, 
    z, m, x,
    hyper, theta, lambda,
    p0, pk)
  {
    gamma <- hyper$gamma
    alpha <- hyper$alpha
    kappa <- hyper$kappa

    gamma_new <- gamma
    hyper$gamma <- gamma_new
    hyper
  }