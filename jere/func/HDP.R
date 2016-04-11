### Note: J is fixed a priori

HDP <- 
  function(
    # data & restaurant number
    x, nJ = 10,
    # DP concentration/base-measure
    gamma = 1, alpha = 1, lik_ji,
    # initialization
    z_init = NULL, m_init = NULL, pi_init = NULL,
    # MCMC parameter
    iter_max = 1000, 
    # utility
    J_init = c("random", "kmeans")[1]
  ){
    #### 0. Global par set up ####
    n <- nrow(x)
    # restaurant mbrshp, either random or kmeans
    J <- # primal rep, n x 1 vector
      switch(J_init, 
             random = sample(1:nJ, n, replace = TRUE),
             kmeans = kmeans(x, nJ)$cluster)
    nJ <- max(J)
    J_idx <- # dual rep, length n list
      lapply(1:nJ, function(i) which(J == i))
    
    #### 1. Initialization ####
    #### > 1.1 default initialization (ground-zero) ====
    if (is.null(z_init)) z_init <- rep(1, n)
    if (is.null(m_init)) m_init <- 0
    if (is.null(pi_init)) pi_init <- c(1, gamma)
    
    #### > 1.2 initialize container ====
    z_iter <- array(NA, dim = c(iter_max, n)) # dish assignment per obs
    m_iter <- vector("list", iter_max) # table count per dish
    pi_iter <- vector("list", iter_max) # dish distribution
    
    z_iter[1, ] <- z_init
    m_iter[[1]] <- m_init
    pi_iter[[1]] <- pi_init
    
    #### 2. Iterate & Visualize ####
    for (ii in 1:iter_max){
      #### 2.1 z: dish assignment per obs ====
      sample_
      
      #### 2.2 m: table count per dish ====
      
      
      #### 2.3 pi: dish distribution ====
      
      
    }
    
    
    #### 3. Return ####
  }