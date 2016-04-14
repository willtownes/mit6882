### Note: J is fixed a priori
HDP <- 
  function(
    # data, restaurant number, dish number
    x, nJ = 1, nK_init = 2,
    # DP concentration/base-measure
    gamma = 1, alpha = NULL, lik_func,
    # initialization
    z_init = NULL, m_init = NULL, pi_init = NULL,
    # MCMC parameter
    iter_max = 1e4,
    # utility
    J_init = c("random", "kmeans")[1]
  ){
    #### 0. Global par set up ####
    n <- nrow(x)
    # restaurant mbrshp, either random or kmeans
    if (is.character(J_init)){
      J <- # primal rep, n x 1 vector
        switch(J_init, 
               random = sample(1:nJ, n, replace = TRUE),
               kmeans = kmeans(x, nJ)$cluster)
    } else J <- J_init
    
    nJ <- max(J)
    J_idx <- # dual rep, length n list
      lapply(1:nJ, function(i) which(J == i))
    
    #### 1. Initialization ####
    #### > 1.1 default initialization (ground-zero) ====
    if (is.null(z_init)) 
      z_init <- 
      sample(1:nK_init, n, replace = TRUE)
    if (is.null(m_init)) 
      m_init <- rep(1, nK_init)
    if (is.null(pi_init)) 
      pi_init <- c(m_init, gamma)/sum(c(m_init, gamma))
    if (is.null(alpha)) 
      alpha <- ceiling(mean(m_init)/2)
    
    #### > 1.2 initialize container ====
    z_iter <- array(NA, dim = c(iter_max, n)) # dish assignment per obs
    m_iter <- vector("list", iter_max) # table count per dish
    pi_iter <- vector("list", iter_max) # dish distribution
    
    #### 2. Iterate & Visualize ####
    z_cur <- z_init
    m_cur <- m_init
    pi_cur <- pi_init
    
    pb <- txtProgressBar(1, iter_max, style = 3)
    for (ii in 1:iter_max){
      setTxtProgressBar(pb, ii)
      #### 2.1 z: dish assignment per obs ====
      z_new <-
        sample_hdp_z(
          x, lik_func, gamma, alpha, J,
          z_cur, m_cur, pi_cur)
      
      #remove empty category
      if(length(setdiff(1:max(z_new), unique(z_new))) > 0){
        cat_idx <- setdiff(1:max(z_new), unique(z_new))
        z_new <- # re-order category
          as.factor(z_new) %>% as.numeric()
        m_cur <- m_cur[-cat_idx]
        pi_cur <- pi_cur[-cat_idx]
      }
      
      #### 2.2 m: table count per dish ====
      m_new <- 
        sample_hdp_m(
          x, lik_func, gamma, alpha, J,
          z_new, m_cur, pi_cur)
      
      #### 2.3 pi: dish distribution ====
      pi_new <- 
        sample_hdp_pi(
          x, lik_func, gamma, alpha, J,
          z_new, m_new, pi_cur)
      
      #### 2.4 update temp par (*_cur) 
      z_iter[ii, ] <- z_cur <- z_new
      m_iter[[ii]] <- m_cur <- m_new
      pi_iter[[ii]] <- pi_cur <- pi_new
      
      if (ii %% 10 == 1){
        plot(x, col = z_cur, 
             main = paste0(ii, ", K = ", length(m_cur))
        )
      }
    }
    
    # discard burn-in 
    iter_burn <- min(round(ii * 0.2), 200)
    z_iter <- z_iter[iter_burn:ii, ]
    m_iter <- m_iter[iter_burn:ii]
    pi_iter <- pi_iter[iter_burn:ii]
    
    #### 3. Return ####
    sim_mat_emp <- 
      lapply(1:1000, 
             function(i)
               outer(z_iter[i, ], z_iter[i, ], "==")
      ) %>% Reduce("+", .) %>% 
      divide_by(nrow(z_iter))
    image(t(sim_mat_emp[nrow(sim_mat_emp):1, ]))
    
  }