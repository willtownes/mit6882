require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func/")
lik_func <- lik_gauss_hdp
data_type <- c("HDP", "manual")[2]

if (data_type == "HDP"){
  set.seed(400)
  dat <- 
    HDP_gauss(# 3 cluster, 1 sub-cluster (per cl), 
      gamma = 10, alpha = 0.5,
      J = 1, K = 10,
      n = 500, d = 2)
} else {
  sigma <- 0.01
  cluster_num <- 5
  each_size <- 50
  
  dat <- vector("list", 2) %>% set_names(c("X", "Z"))
  dat$X <- 
    expand.grid(1:cluster_num, 1:cluster_num) %>% 
    t %>% as.data.frame %>% 
    lapply(function(mean) 
      rmvnorm(each_size, mean, sigma*diag(2))) 
  dat$Z <- rep(1:cluster_num^2, each = each_size)
} 

x <- dat %>% extract2(1) %>% do.call(rbind, .)
z <- dat[[2]] %>% as.factor %>% as.numeric
x <- x[order(z), ]
plot(x)

HDP(
  # data, restaurant number, dish number
  x, nJ = 1, nK_init = 5,
  # DP concentration/base-measure
  gamma = 0.1, alpha = 0.5, lik_func,
  # initialization
  z_init = NULL, m_init = NULL, pi_init = NULL,
  # MCMC parameter
  iter_max = 1e3,
  # utility
  J_init = c("random", "kmeans")[1]
)
