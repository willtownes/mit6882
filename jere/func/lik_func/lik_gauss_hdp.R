library(mvtnorm)
# likelihood for predictive dist x_{k}^{-(ji)} | x_{k}^{-(ji)} 
# from gaussian likelihood: 
#     x ~ N(u, I)
#     u ~ N(0, sigma * I)

lik_gauss_hdp <- 
  function(x, z, m,
           # misc data information
           n, d, K,
           # prior par 
           sigma_dat = 0.01, 
           sigma = NULL, log = FALSE){
    # input:
    # > x (n x p):      data
    # > z (1 x n):      cluster assignment 
    # > sigma (1 x 1):  prior variance
    # output:
    # > lik [n x (K+1)]: predictive prob for each x^{ji}
    # > n [n x (K+1)]: sample size for each x^{ji}
    if (is.null(sigma)) sigma <- 100 #K/n
    # prevent missing category scenario
    K_cat <- sort(unique(z))
    
    # 1. calculate mean & variance weights 
    x_k <- # obs assigned to each cluster
      lapply(K_cat, 
             function(z_cat) 
               # prevent 1 row scenario
               matrix(x[which(z == z_cat), ], ncol = d)
      )
    
    # n_k: num obs per cluster
    n_k <- sapply(x_k, nrow)
    n_ji <- lapply(
      1:n, function(i){
        c(n_k - (z[i] == K_cat), 0)
      }) %>% do.call(rbind, .)
    
    # mean/var weights, n obs with (K+1) x 2 list
    mult_ji <- 
      lapply(
        1:n, function(i){
          n_obs <- n_ji[i, ]
          list(mean = sigma * n_obs / (sigma * n_obs + sigma_dat), 
               var =  sigma_dat * (1 + sigma / (sigma * n_obs + sigma_dat))
          )
        })
    
    # 2. calculate mean_ji
    sum_k <- lapply(x_k, colSums) %>% do.call(rbind, .)
    mean_ji <- # mean for ji^th obs
      lapply(1:n, 
             function(i) # (K+1) x d matrix
               if (n_k[[ z[i] ]] == 1){
                 # deal with the n=1 case
                 (sum_k - 
                    matrix((z[i] == K_cat)) %*% t(x[i, ]) ) %>%
                   # add empty "mean" for new class
                   rbind(t(rep(0, d)))
               } else {
                 ((#t(sum_k[[ z[i] ]])[rep(1, K), ] - 
                   sum_k - matrix((z[i] == K_cat)) %*% t(x[i, ]))/
                    (n_k - (z[i] == K_cat))) %>%
                   # add empty "mean" for new class
                   rbind(t(rep(0, d)))
               }
      )
    
    # check if mean are calculated correctly
    # mean_K <-
    #   lapply(1:K,
    #          function(k){
    #            lapply(mean_ji, function(x) x[k, ]) %>%
    #              do.call(rbind, .)
    #          }
    #   )
    # plot(x, col = "grey")
    # text(mean_K[[1]], label = "1")
    # text(mean_K[[2]], label = "2")
    # text(mean_K[[3]], label = "3")
    
    # 3. calculate likelihood
    dmvnorm_V <- 
      Vectorize(
        function(i, k)
          dmvnorm(x[i, ],
                  mult_ji[[i]]$mean[k] * mean_ji[[i]][k, ], 
                  diag(rep(mult_ji[[i]]$var[k], d)))
      )
    lik <- 
      outer(X = 1:n, Y = 1:(K+1), dmvnorm_V)
    if (log) lik <- log(lik)
    
    # # examine likelihood allocation
    # k_id <- 3
    # plot(x, pch = 19, 
    #      col = rgb(0,0,0, lik[, k_id]/max(lik[, k_id])))
    
    # return
    list(lik = lik, n_ji = n_ji)
  }