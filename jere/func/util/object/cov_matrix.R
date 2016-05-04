## Generate covariance matrix based on given cor and cov parameter
# actual cor and var are generated using rexp

cov_matrix <- function(d, cor_mean, var_mean){
  # initialize matrix with random cor value
  covMat <- 
    rexp(d^2, 1/cor_mean) %>% matrix(ncol = d)
  
  # fill in variance
  diag(covMat) <- rexp(d, 1/var_mean)
  
  # symmetrify (take lower.tri)
  covMat[lower.tri(covMat)] <- t(covMat)[lower.tri(covMat)]
  
  # make matrix PD through ridge correction
  correct_val <- 
    eigen(covMat, only.values = TRUE)$values %>% 
    min %>% min(0, .)
  diag(covMat) <- diag(covMat) - correct_val + 1e-5
  
  # return
  covMat
}