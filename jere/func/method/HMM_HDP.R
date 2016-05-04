### Note: Sticky HDP applied to HMM setting
HMM_HDP <- 
  function(
    # data, restaurant number, dish number
    x, K = 10,
    # base measure & parameters
    log_lik_func = NULL, lambda_init = NULL, 
    sample_lambda = FALSE, lambda_sampler = NULL, 
    # DP concentration/base-measure
    gamma = 1, alpha = NULL, kappa = NULL, 
    # initialization
    z_init = NULL, m_init = NULL, 
    p0_init = NULL, pk_init = NULL,
    # MCMC parameter
    iter_max = 2e3, iter_update = 10, verbose = FALSE)
  {
    # K:           max number of states
    # log_lik_func:    log emission dist likelihood
    # lambda_init: initial value of emission parameters, 1:K list
    #
    # iter_update: num of iteration to wait before update plot
    
    #### 0. Global par set up ####
    n <- nrow(x)
    d <- ncol(x)
    
    if (is.null(theta_init)|is.null(log_lik_func)) 
      stop("emission parameters 'theta_init' and 'log_lik_func' must be supplied")
    if (sample_lambda){
      if (is.null(lambda_sampler))
        stop("sample_lambda = TRUE but 'lambda_sampler' not supplied")
      if (length(lambda_sampler) != K)
        stop("'lambda_sampler' should be a list of length K")
    }

    #### 1. Initialization ####
    #### > 1.1 default initialization (random assignment) ====
    if (is.null(z_init)) # random assign 
      z_init <- sample(1:K, n, replace = TRUE)
    if (is.null(m_init)) # equal table number
      m_init <- rep(1, K)
    if (is.null(p0_init)) # equal global dish prob
      p0_init <- rep(1/K, K)
    if (is.null(pk_init)) # equal restaurant-dish prob
      pk_init <- array(1/K, dim = c(K, K))
    
    if (is.null(alpha)) alpha <- 1
    if (is.null(kappa)) kappa <- 0.5
    
    #### > 1.2 initialize container ====
    z_iter <- array(NaN, dim = c(iter_max, n)) # state assignment per obs
    m_iter <- array(NaN, dim = c(iter_max, K)) # table count per dish
    p0_iter <- array(NaN, dim = c(iter_max, K)) # global state distribution
    pk_iter <- array(NaN, dim = c(iter_max, K, K)) # state transition distribution
    lambda_iter <- vector("list", length = iter_max)
    
    # optionally, sample hyperparameters
    gamma_iter <- array(NaN, dim = c(iter_max, 1))
    alpha_iter <- array(NaN, dim = c(iter_max, 1))
    kappa_iter <- array(NaN, dim = c(iter_max, 1))
    
    #### 2. Iterate & Visualize ####
    z_cur <- z_init
    m_cur <- m_init
    p0_cur <- p0_init
    pk_cur <- pk_init
    
    gamma_cur <- gamma
    alpha_cur <- alpha
    kappa_cur <- kappa
    lambda_cur <- lambda_init
    
    pb <- txtProgressBar(1, iter_max, style = 3)
    for (ii in 1:iter_max){
      setTxtProgressBar(pb, ii)
      #### 2.1 z: state assignment per obs ====
      z_new <-
        sample_hmm_z(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_cur, 
          z = z_cur, m = m_cur, 
          p0 = p0_cur, pk = pk_cur,
          verbose)

      #### 2.2 m: table count per state ====
      m_list <-
        sample_hmm_m(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_cur, 
          z = z_new, m = m_cur, 
          p0 = p0_cur, pk = pk_cur)
      
      m_new <- m_list$m
      n_jk <- m_list$n
      
      #### 2.3 p0 & pk: state-transition ====
      p0_new <-
        sample_hmm_p0(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_cur, 
          z = z_new, m = m_new,
          p0 = p0_cur, pk = pk_cur, n_jk = n_jk)
 
      pk_new <-
        sample_hmm_pk(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_cur, 
          z = z_new, m = m_new,
          p0 = p0_new, pk = pk_cur, n_jk = n_jk)     
      
      #### 2.4 theta: emission par ====
      # to be plugged in
      lambda_new <- 
        sample_hmm_lambda(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_cur, 
          z = z_new, m = m_new,
          p0 = p0_new, pk = pk_new, n_jk = n_jk)       
      
      #### 2.5 hyperparameters ====
      gamma_new <- # to do 
        sample_hmm_gamma(
          x, log_lik_func, 
          gamma = gamma_cur, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_new, 
          z = z_new, m = m_new,
          p0 = p0_new, pk = pk_new, n_jk = n_jk)          
      
      alpha_new <- # to do
        sample_hmm_alpha(
          x, log_lik_func, 
          gamma = gamma_new, alpha = alpha_cur, 
          kappa = kappa_cur, lambda = lambda_new, 
          z = z_new, m = m_new,
          p0 = p0_new, pk = pk_new, n_jk = n_jk)          
      
      kappa_new <- # to do
        sample_hmm_kappa(
          x, log_lik_func, 
          gamma = gamma_new, alpha = alpha_new, 
          kappa = kappa_cur, lambda = lambda_new, 
          z = z_new, m = m_new,
          p0 = p0_new, pk = pk_new, n_jk = n_jk)

      #### 2.z update temp par (*_cur) 
      z_iter[ii, ] <- z_cur <- z_new
      m_iter[ii, ] <- m_cur <- m_new
      p0_iter[ii, ] <- p0_cur <- p0_new
      pk_iter[ii, ,] <- pk_cur <- pk_new
      
      gamma_iter[ii, ] <- gamma_cur <- gamma_new
      alpha_iter[ii, ] <- alpha_cur <- alpha_new
      kappa_iter[ii, ] <- kappa_cur <- kappa_new
      
      
      if (ii %% iter_update == 1){
        plot(x[, 1], col = z_cur, 
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
    z_iter[length(z_iter)]
  }