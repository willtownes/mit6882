library(MASS)
datgen_2d <- 
  function(n, d = 2){
  # n: number of clusters
  mu <- cbind(runif(n), runif(n))
  sig <- 1e-3
  N <- rpois(n, 40)
  
  sample <- 
    lapply(1:n, 
           function(i) 
             mvrnorm(N[i], mu[i,], diag(rep(sig, d)))
    )
  do.call(rbind, sample)
}