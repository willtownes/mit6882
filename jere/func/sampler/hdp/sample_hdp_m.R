sample_hdp_m <- 
  function(
    # data & likelihood
    x, lik_func, 
    # CRF parameter, prior
    gamma, alpha, J, 
    # CRF parameter, data
    z, m, pi){
    # sample_hdp_m: sample table count for dish by simulating t_ji
    ######################################################
    # sample z, the customer-specific dish assignment
    # input:
    # > lik_func: likelihood function return lik and n_ji
    # > z:  dish assignment for obs
    # > m:  table count per parameter
    # output:
    # > m: updated table count
    ######################################################    
    
    nJ <- max(J)
    K_cat <- sort(unique(z))
    K <- length(K_cat)
    
    #### 1. sample t_ji ====
    # 1.1 initialize data, customer count per dish
    x_j <- 
      lapply(1:nJ, function(j) x[J==j, ])
    n_jk <- 
      # obs count per dish per restaurant
      tapply(z, J,
             function(z_j){
               sapply(K_cat, 
                      function(k) 
                        sum(z_j==k)) %>% set_names(K_cat)
             })
    
    # 1.2 compute m_jk
    m_jk <- 
      lapply(1:nJ, 
             # for each restaurant
             function(j)
             { # compute per-restaurant, dish-spec table count
               sapply(1:K, 
                      function(k){
                        n <- n_jk[[j]][k]
                        # vector of prob, increment n [(1:n) -1], then sample
                        ((alpha * pi[k])/(0:n + alpha * pi[k])) %>%
                          # sample new t_jk
                          sapply(function(p) rbinom(1, 1, p)) %>%
                          # m_jk
                          sum
                      }
               )
             }
      )
    
    # 1.3 compute m_k, return
    m <- do.call(cbind, m_jk) %>% rowSums()
  }