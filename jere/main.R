source("./func/util/source_Dir.R")
sourceDir("./func/")

set.seed(100)
x <- 
  gauss_2d(# 3 cluster, 1 sub-cluster (per cl), 
    k = 3, t = 1, d = 2, 
    sig_k = 1e-3, sig_t = 5e-4, 
    mean_obs = 100) %>% scale

plot(x)



lik_func <- lik_gauss_mean
