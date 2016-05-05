### Note: Sticky HDP applied to HMM setting
HMM_HDP <- 
  function(
    # data, restaurant number, dish number
    y, K = 2,
    # emission measure & parameters
    log_lik_func = NULL, 
    sample_emiss = FALSE, 
    theta_init = NULL, theta_sampler = NULL, 
    lambda_init = NULL, lambda_sampler = NULL,
    x_init = NULL, x_sampler = NULL,
    # DP concentration/base-measure
    sample_hyper = FALSE,
    gamma = 1, alpha = NULL, kappa = NULL, 
    # initialization
    z_init = NULL, m_init = NULL, 
    p0_init = NULL, pk_init = NULL,
    # MCMC parameter
    iter_max = 1e3, iter_update = 10, verbose = FALSE)
  {
    # INPUT:
    # K:              max number of states
    # log_lik_func:   log emission dist likelihood
    # x_init:         initial for latent obs including time 0
    #                 (T+1) x d df
    # theta_init:     initial for state-specific emission parameters, 
    #                 1:K list
    # lambda_init:    initial for state-invariant emission par, list
    # iter_update:    num of iteration to wait before update plot
    
    #### 0. Global par set up ####
    T <- nrow(y)
    d <- ncol(y)
    
    if (is.null(x_init)|is.null(theta_init)|is.null(lambda_init)) 
      stop("initial values for x/theta/lambda must be supplied")
    if (is.null(theta_sampler)|is.null(lambda_sampler)) 
      stop("emission samplers for theta/lambda must be supplied")
    if (length(theta_init) != K)
      stop("'theta_sampler' should be a list of length K")
    if (sample_emiss){
      if (is.null(theta_sampler)|is.null(lambda_sampler))
        stop("sample_emiss = TRUE but lambda/theta sampler not supplied")
    }
    
    #### 1. Initialization ####
    #### > 1.1 default initialization (random assignment) ====
    if (is.null(z_init)) # random assign 
      z_init <- sample(1:K, T, replace = TRUE)
    if (is.null(m_init)) # equal table number
      m_init <- rep(1, K)
    if (is.null(p0_init)) # equal global dish prob
      p0_init <- rep(1/K, K)
    if (is.null(pk_init)) # equal restaurant-dish prob
      pk_init <- array(1/K, dim = c(K, K))
    
    if (is.null(alpha)) alpha <- 1
    if (is.null(kappa)) kappa <- 0.5
    
    #### > 1.2 initialize container ====
    z_iter <- array(NaN, dim = c(iter_max, T)) # state assignment per obs
    m_iter <- array(NaN, dim = c(iter_max, K)) # table count per dish
    p0_iter <- array(NaN, dim = c(iter_max, K)) # global state distribution
    pk_iter <- array(NaN, dim = c(iter_max, K, K)) # state transition distribution
    
    # optionally, sample emission parameters
    if (sample_emiss){
      theta_iter <- vector("list", length = iter_max)
      lambda_iter <- vector("list", length = iter_max)
    }
    # optionally, sample hyperparameters
    if (sample_hyper){
      gamma_iter <- array(NaN, dim = c(iter_max, 1))
      alpha_iter <- array(NaN, dim = c(iter_max, 1))
      kappa_iter <- array(NaN, dim = c(iter_max, 1))
    }
    
    #### 2. Iterate & Visualize ####
    x_cur <- x_init
    z_cur <- z_init
    m_cur <- m_init
    p0_cur <- p0_init
    pk_cur <- pk_init
    
    hyper_cur <- 
      list(gamma = gamma, alpha = alpha, kappa = kappa)
    theta_cur <- theta_init
    lambda_cur <- lambda_init
    
    pb <- txtProgressBar(1, iter_max, style = 3)
    for (ii in 1:iter_max){
      setTxtProgressBar(pb, ii)
      
      #### 2.1 x: pseudo obs for SLDS ====
      x_new <-
        x_sampler(
          y, log_lik_func, 
          z = z_cur, m = m_cur, x = x_cur,
          p0 = p0_cur, pk = pk_cur,          
          theta = theta_cur, lambda = lambda_cur,
          hyper = hyper_cur
        )
      
      #### 2.1 z: state assignment per obs ====
      z_new <-
        sample_hmm_z(
          y, log_lik_func, 
          z = z_cur, m = m_cur, x = x_new,
          p0 = p0_cur, pk = pk_cur,          
          theta = theta_cur, lambda = lambda_cur,
          hyper = hyper_cur, 
          verbose = verbose
        )
      
      #### 2.2 m: table count per state ====
      m_list <-
        sample_hmm_m(
          y, log_lik_func, 
          z = z_new, m = m_cur, x = x_new,
          p0 = p0_cur, pk = pk_cur,          
          theta = theta_cur, lambda = lambda_cur,
          hyper = hyper_cur
        )
      
      m_new <- m_list$m
      hyper_cur$n_jk <- m_list$n
      
      #### 2.3 p0 & pk: state-transition ====
      p0_new <-
        sample_hmm_p0(
          y, log_lik_func, 
          z = z_new, m = m_new, x = x_new,
          p0 = p0_cur, pk = pk_cur,          
          theta = theta_cur, lambda = lambda_cur,
          hyper = hyper_cur
        )
      
      pk_new <-
        sample_hmm_pk(
          y, log_lik_func, 
          z = z_new, m = m_new, x = x_new,
          p0 = p0_new, pk = pk_cur,          
          theta = theta_cur, lambda = lambda_cur,
          hyper = hyper_cur
        )     
      
      #### 2.4 theta & lambda: emission par ====
      if (sample_emiss){
        # to be plugged in
        theta_new <- 
          theta_sampler(
            y, log_lik_func, 
            z = z_new, m = m_new, x = x_new,
            p0 = p0_new, pk = pk_new,          
            theta = theta_cur, lambda = lambda_cur,
            hyper = hyper_cur)       
        
        lambda_new <- 
          lambda_sampler(
            y, log_lik_func, 
            z = z_new, m = m_new, x = x_new,
            p0 = p0_new, pk = pk_new,          
            theta = theta_new, lambda = lambda_cur,
            hyper = hyper_cur) 
      }
      
      #### 2.5 hyperparameters ====
      # sample new hyperpar & insert into the list 'hyper'
      if (sample_hyper){
        hyper_new <- # to do
          sample_hmm_gamma(
            y, log_lik_func, 
            z = z_new, m = m_new, x = x_new,
            p0 = p0_new, pk = pk_new,          
            theta = theta_new, lambda = lambda_new,
            hyper = hyper_cur)          
        
        hyper_new <- # to do
          sample_hmm_alpha(
            y, log_lik_func, 
            z = z_new, m = m_new, x = x_new,
            p0 = p0_new, pk = pk_new,          
            theta = theta_new, lambda = lambda_new,
            hyper = hyper_new)          
        
        hyper_new <- # to do
          sample_hmm_kappa(
            y, log_lik_func, 
            z = z_new, m = m_new, x = x_new,
            p0 = p0_new, pk = pk_new,          
            theta = theta_new, lambda = lambda_new,
            hyper = hyper_new)
      }
      
      #### 2.z update temp par (*_cur) 
      z_iter[ii, ] <- z_cur <- z_new
      m_iter[ii, ] <- m_cur <- m_new
      p0_iter[ii, ] <- p0_cur <- p0_new
      pk_iter[ii, ,] <- pk_cur <- pk_new
      
      if (sample_emiss){
        theta_iter[[ii]] <- theta_cur <- theta_new
        lambda_iter[[ii]] <- lambda_cur <- lambda_new
      }
      
      if (sample_hyper){
        hyper_cur <- hyper_new
        gamma_iter[ii, ] <- hyper_cur$gamma
        alpha_iter[ii, ] <- hyper_cur$alpha
        kappa_iter[ii, ] <- hyper_cur$kappa
      }
      
      if (ii %% iter_update == 0){
        plot(y, col = z_cur, pch = 19,
             main = paste0(ii, ", K = ", length(m_cur))
        )
      }
    }
    
    # discard burn-in 
    iter_burn <- min(round(ii * 0.2), 200)
    z_iter <- z_iter[iter_burn:ii, ]
    m_iter <- m_iter[iter_burn:ii, ]
    p0_iter <- p0_iter[iter_burn:ii, ]
    pk_iter <- pk_iter[iter_burn:ii, , ]
    
    #### 3. Return ####
    z_iter
  }