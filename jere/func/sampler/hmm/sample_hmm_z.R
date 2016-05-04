sample_hmm_z <- 
  function(
    # data & likelihood
    x, log_lik_func, 
    # CRF parameter, prior
    gamma, alpha, kappa, lambda,
    # CRF parameter, data
    z, m, p0, pk, 
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
    n <- nrow(x)
    d <- ncol(x)
    K <- length(m)

    #### 1. compute backward message ####
    # a state*time message matrix, 
    # (last col reserved for initial value)
    # SLOOOOWWWWWWWWW
    m_KT <- 
      backward_message(x, lambda, pk, 
                       log_lik_func, 
                       verbose)
    
    #### 2. assemble forward likelihood  ####
    f_KT <- # a state*time prob matrix
      forward_prob(x, z, lambda, pk, 
                   log_lik_func, m_KT, 
                   verbose)
    
    #### 3. perform forward sampling  ####
    z_new <- 
      apply(f_KT, 2, function(prob) rMulti(1, prob))
    
    #### 4. clean up n return ####
    return(z_new)
  }