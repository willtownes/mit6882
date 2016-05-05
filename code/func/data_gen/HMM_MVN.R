library(MASS)
library(magrittr)

HMM_MVN <- 
  function(Z_T, 
           d_gauss = 1, cor_gauss = 0.2, var_gauss = 1)
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
    
    # 2. generate data
    X_T <- 
      lapply(Z_T, 
             function(z) 
               rmvnorm(1, theta_0[[z]], Sigma_0[[z]])
      ) %>% do.call(rbind, .)
    
    # return sample
    list(X = X_T, Z = Z_T, 
         Theta = 
           lapply(1:K,
                  function(i) 
                    list(mean = theta_0[[i]], cov = Sigma_0[[i]])
           )
    )
  }