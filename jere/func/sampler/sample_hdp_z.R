sample_hdp_z <- 
  function(
    # data & likelihood
    x, lik_func, 
    # CRF parameter, prior
    gamma, alpha, J, 
    # CRF parameter, data
    z, m, pi)
    ######################################################
    # sample z, the customer-specific dish assignment
    # input:
    # > lik_func: likelihood function return lik and n_ji
    # > z:  dish assignment for obs
    # > m:  table count per parameter
    # output:
    # > z: updated z assigment
    ######################################################    
  {
    # 0. initialize parameters
    n <- nrow(x)
    d <- ncol(x)
    K <- length(m)
    
    # 1. calculate likelihood
    # invoke custom log-lik func, SLOW
    lik_out <- 
      lik_func(x, z, m, # data & parameter
               n, d, K) # misc information
    lik_val <- lik_out$lik # log likhd
    n_ji <- lik_out$n_ji
    
    # 2. calculate weights
    weight <- n_ji + t((alpha * pi))[rep(1, n), ]
    
    # 3. assemble likelihood
    z_prob <-
      (lik_val * weight) %>%
      apply(1, function(row) 
        log(row) - log(sum(row))) %>% 
      exp %>% t
    
    # 4. sample
    z_new <- 
      apply(z_prob, 1,
            function(prob)
              which(rmultinom(1, 1, prob)>0)
      )
    
    return(z_new)
  }