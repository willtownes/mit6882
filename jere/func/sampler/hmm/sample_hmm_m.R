library(plyr)

sample_hmm_m <- 
  function(
    y, log_lik_func, 
    z, m, x,
    p0, pk,
    theta, lambda,
    hyper
  )
  {
    # sample_hmm_m: sample table count for dish by simulating t_ji
    ######################################################
    # input:
    # > z:  dish assignment for obs
    # > m:  table count per parameter
    # output:
    # > m_new: updated table count
    ######################################################    
    
    #### 0. initialize parameters ####
    n <- nrow(y)
    d <- ncol(y)
    K <- length(m)
    
    gamma <- hyper$gamma
    alpha <- hyper$alpha
    kappa <- hyper$kappa
    
    #
    m_jk <- matrix(NaN, nrow = K, ncol = K)
    n_jk <- 
      # long format
      cbind(z[-n], z[-1]) %>% as.data.frame %>%
      plyr::count(vars = c("V1", "V2")) %>%
      # matrix format
      (function(tb){
        # create empty matrix
        out <- matrix(0, nrow = K, ncol = K)
        # fill in entries from long format
        for (i in 1:nrow(tb)){
          out[tb[i, 1], tb[i, 2]] <- tb[i, 3]
        }
        out
      })
    
    
    #### 1. sample naive table count m' ####
    p_mat <- 
      matrix(rep(p0, each = K), nrow = K) * alpha +
      diag(K) * kappa
    
    m_jk <- 
      (function(j, k){
        m_jk <- 0
        for (n_incr in 1:n_jk[j, k]){
          p <- p_mat[j, k]/(n_incr + p_mat[j, k])
          m_jk <- m_jk + rbinom(1, 1, p)
        }
        m_jk
      }) %>% Vectorize %>% outer(1:K, 1:K, .)
    
    #### 2. sample override variable w ####
    p_w <- kappa/(kappa + alpha * p0)
    w <- sapply(
      1:K, function(k) rbinom(1, m_jk[k, k], p_w[k])) 
    
    #### 3. sample effective table count m ####
    m_jk <- m_jk - diag(w)
    
    #### 4. assemble n return ####
    return(list(m = colSums(m_jk), n = n_jk))
  }