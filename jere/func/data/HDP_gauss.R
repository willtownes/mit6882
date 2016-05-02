library(MASS)
library(magrittr)

HDP_gauss <- 
  function(
    gamma = 100, alpha = 0.5,
    J = 1, K = 10,
    n = 100, d = 2
  )
  {
    # K: number of dishes
    # J: number of restaurants
    # n: number of customers per restaurant
    
    # 1 Draw G0, a dataframe of theta - beta
    theta_0 <- matrix(runif(K*d, -1, 1), ncol = d)
    beta_0 <- rGEM(K, alpha = gamma)
    
    # 2 Draw GJ, draw K*10 many times
    G_J <- 
      (function(j){
        idx_j <- rMulti(K*10, beta_0)
        beta_j <- 
          rGEM(K*10, alpha = alpha) %>% 
          tapply(idx_j, sum)
        beta_j
      }) %>% lapply(1:J, .)
    
    # visualize distribution
    # lapply(G_J, function(x) plot(names(x), x))
    
    # 3. generate z_{ji} from Gj
    Z_J <- 
      lapply(G_J, function(beta)
        names(beta)[rMulti(n, beta)] %>% 
          as.numeric
      )
    
    # 4. generate data
    X_J <- 
      lapply(Z_J, function(z_j){
        sapply(z_j, function(z) 
          rmvnorm(1, theta_0[z, ], 0.01*diag(rep(1, d)))
        ) %>% t
      })
    
    # plot(do.call(rbind, X_J), col = unlist(Z_J))
    # points(theta_0, pch = 15, col = "white", cex = 4)
    # text(theta_0, labels = 1:10, col = 2, cex = 2)
    
    # return sample
    list(X_J, z = unlist(Z_J))
  }