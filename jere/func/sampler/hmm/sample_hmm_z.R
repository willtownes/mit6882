sample_hmm_z <- 
  function(
    y, log_lik_func, 
    z, m, x,
    p0, pk,
    theta, lambda,
    hyper,
    verbose = FALSE
    )
    ######################################################
    # sample z, the obs-spec state assignment in hmm 
    #           using forward-backward msg passing
    # input:
    # > log_lik_func: likelihood function return lik and n_ji
    # > z:  dish assignment for obs
    # > m:  table count per parameter
    # output:
    # > z: updated z assigment
    ######################################################    
  {
    #### 0. initialize parameters ####
    T <- nrow(y)
    d <- ncol(y)
    K <- length(m)

    #### 1. compute backward message ####
    # a state*time message matrix, 
    # (last col reserved for initial value)
    # SLOOOOWWWWWWWWW
    m_KT <- 
      backward_message(y, x, z,
                       theta, pk, log_lik_func, verbose)
    
    #### 2. assemble forward likelihood  ####
    f_KT <- # a state*time prob matrix
      forward_prob(y, x, z, m_KT, 
                   theta, pk, log_lik_func, verbose)
    
    #### 3. perform forward sampling  ####
    z_new <- 
      apply(f_KT, 2, function(prob) rMulti(1, prob))
    
    #### 4. clean up n return ####
    return(z_new)
  }