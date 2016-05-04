library(MASS)
library(magrittr)

HMM_gauss <- 
  function(
    gamma = 100, alpha = 10, kappa = 1, 
    K = 100, T = 100, 
    d_gauss = 1, cor_gauss = 0.2, var_gauss = 1
  )
  {
    # Hidden Markov Model with MV Gaussian Emissions
    
    # K: number of states
    # T: number of time points
    # d: dimension of MVN emission
    
    # 0 candidate parameters
    theta_0 <- 
      lapply(1:K, 
             function(k) rep(k, d_gauss))
    Sigma_0 <- 
      lapply(1:K, 
             function(k) 
               cov_matrix(d = d_gauss, 
                          cor_mean = cor_gauss, 
                          var_mean = var_gauss)
      )
    
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
    
    # 4. generate data
    X_T <- 
      lapply(Z_T, 
             function(z) 
               rmvnorm(1, theta_0[[z]], Sigma_0[[z]])
      ) %>% do.call(rbind, .)
    
    # plot(X_T)
    
    # return sample
    list(X = X_T, Z = Z_T, 
         Theta = 
           lapply(1:K,
                  function(i) 
                    list(mean = theta_0[[i]], cov = Sigma_0[[i]])
           )
    )
  }