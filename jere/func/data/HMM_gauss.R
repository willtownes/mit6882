library(MASS)
library(magrittr)

HMM_gauss <- 
  function(
    gamma = 100, alpha = 10, kappa = 1, 
    K = 100, T = 100
  )
  {
    # Hidden Markov Model with Gaussian Emissions
    
    # K: number of states
    # t: number of time points
    
    # 1 Draw G_0, a dataframe of theta - beta
    theta_0 <- runif(K)
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
    
    # 4. generate data
    X_T <- 
      sapply(Z_T, 
             function(z) rnorm(1, theta_0[z], 0.01))
    
    # plot(1:T, X_T)

    # return sample
    list(X = X_T, Z = Z_T)
  }