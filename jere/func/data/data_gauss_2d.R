library(MASS)
library(magrittr)

gauss_2d <- 
  function(k = 10, t = 3, d = 2, 
           sig_k = 1e-2, sig_t = 1e-3)
  {
    # k: number of clusters
    # n: number of sub-clusters in each cluster
    
    # 1. determine cluster/subcluster mean
    mu_k <- cbind(runif(k), runif(k))
    mu_t <- # generate per cluster mean
      lapply(1:k, 
             function(i) 
               mvrnorm(t, mu_k[i,], diag(rep(sig_k, d)))
      ) %>% do.call(rbind, .)
    
    # 2. generate data for each sub-cluster
    n_t <- nrow(mu_t)
    N <- rpois(n_t, 40)
    sample <-       
      lapply(1:n_t,
             function(i) 
               mvrnorm(N[i], mu_t[i,], diag(rep(sig_t, d)))
      ) %>% do.call(rbind, .)
  
    # return sample
    sample
  }