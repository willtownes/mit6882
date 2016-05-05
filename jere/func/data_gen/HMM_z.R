library(MASS)
library(magrittr)

HMM_z <- 
  function(
    gamma = 100, alpha = 10, kappa = 1, 
    K = 100, T = 100
  )
  {
    # generate mode allocation for Hidden Markov Model
    
    # K: number of states
    # T: number of time points
    # d: dimension of MVN emission
    
    # 1 Draw G_0, a dataframe of theta - beta
    beta_0 <- rGEM(K, alpha = gamma)
    
    # 2 Draw G_K, draw K*10 many times
    G_K <- 
      (function(j){
        idx_k <- rMulti(K*10, beta_0)
        beta_k <- 
          rGEM(K*10, alpha = alpha) %>% 
          tapply(idx_k, sum)
        # reorder category
        beta_out <- rep(0, K)
        beta_out[as.numeric(names(beta_k))] <- beta_k
        names(beta_out) <- 1:K
        beta_out
      }) %>% lapply(1:K, .)
    
    # visualize distribution
    # lapply(G_K, function(x) plot(names(x), x))
    
    # 3. generate z_{ji} from Gj
    Z_T <- rep(0, T)
    Z_T[1] <- rMulti(1, beta_0)
    for(t in 2:T){
      z_prev <- Z_T[t-1]
      # sticky transition
      beta_prev <- G_K[[ z_prev ]] 
      idx <- G_K
      beta_prev[z_prev] <-  beta_prev[z_prev] + kappa
      beta_prev <- beta_prev/sum(beta_prev)
      #
      Z_T[t] <- rMulti(1, beta_prev)
    }
    
    Z_T
  }