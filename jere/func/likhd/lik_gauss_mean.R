# likelihood for predictive dist x_{k}^{-(ji)} | x_{k}^{-(ji)} 
# from gaussian likelihood: 
#     x ~ N(u, I)
#     u ~ N(0, sigma * I)
library(mvtnorm)

lik_gauss_mean <- 
  function(x, n, K, z,
           # prior par 
           sigma = 1e2){
    # input:
    # > x (n x p):      data
    # > z (1 x n):      cluster assignment 
    # > sigma (1 x 1):  prior variance
    # output:
    # > lik [n x (K+1)]: predictive prob for each x^{ji}
    
    d <- ncol(x)
    # 1. calculate mean & variance weights 
    x_k <- #obs assigned to each cluster
      lapply(1:K, 
             function(z_cat) x[which(z == z_cat), ])
    # n_k: num obs per cluster
    n_k <- lapply(x_z, nrow)
    # mean/var weights, n obs with (K+1) x 2 list
    mult_ji <- 
      lapply(
        1:n, function(i){
          n_ji <- 
            c(n_k[[ z[i] ]] - (z[i] == 1:K), 0)
          list(mean = sigma * n_ji / (sigma * n_ji + 1), 
               var =  1 + sigma / (sigma * n_ji + 1)
          )
        })
    
    
    # 2. calculate mean_ji
    sum_k <- lapply(x_z, colSums)
    mean_ji <- # mean for ji^th obs
      lapply(1:n, 
             function(i) # (K+1) x d matrix
               ((t(sum_k[[ z[i] ]])[rep(1, K), ] - 
                   matrix((z[i] == 1:K)) %*% t(x[i, ]) )/
                  (n_k[[ z[i] ]] - (z[i] == 1:K))) %>%
               # add empty "mean" for new class
               rbind(t(rep(0, d)))
      )
    
    # 3. calculate likelihood
    dmvnorm_V <- Vectorize(
      function(i, k)
        dmvnorm(x[i, ],
                mult_ji[[i]]$mean[k] * mean_ji[[i]][k, ], 
                diag(rep(mult_ji[[i]]$var[k], d)))
    )
    lik <- 
      outer(X = 1:n, Y = 1:(K+1), dmvnorm_V)
    
    # return
    lik
  }