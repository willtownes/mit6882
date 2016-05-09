library(MASS)
library(magrittr)

HMM_SLDS <- function(Z_T, theta, lambda = NULL, 
                     x_0 = NULL,
                     d_y = 2)
{
  # Hidden Markov Model with MV Gaussian Emissions
  
  # Z_T: mode assignment for observations
  # theta: input SLDS parameters, mode specific
  # lambda: input SLDS parameters, not mode specific
  # d_y: dimension of SLDS emission, 
  #     used only if lambda don't contain "C" 
  
  # 0 candidate parameters
  d <- length(theta[[1]]$B)
  T <- length(Z_T)
  
  # 1. generate X
  X_T <- matrix(NaN, nrow = T+1, ncol = d)
  if (is.null(x_0)){
    X_T[1, ] <- 0 # initiate x_0
  } else {
    X_T[1, ] <- x_0
  }
  
  for (t in (1:T) + 1){
    z_t <- Z_T[t-1]
    mu_t <- 
      theta[[z_t]]$A %*% X_T[t-1, ] + theta[[z_t]]$B
    X_T[t, ] <- 
      rmvnorm(1, mu_t, theta[[z_t]]$Sigma)
  }

  # 2. generate Y
  # 2.1 create C and R matrix
  if (!("C" %in% names(lambda))){
    lambda$C <- 
      cbind(diag(d_y), 
            matrix(0, nrow = d_y, ncol = d - d_y)
      )
  }
  if (!("R" %in% names(lambda))){
    lambda$R <- diag(d_y)/100
  }
  
  Y_T <- matrix(NaN, nrow = T, ncol = nrow(lambda$C))
  Y_T <- X_T[-1, ] %*% t(lambda$C) + rmvnorm(T, rep(0, 2), lambda$R)
  
  
  plot(X_T[, 1:2], type = "l")
  points(Y_T)
  
  # return sample
  list(Y = Y_T, X = X_T, Z = Z_T, Theta = theta)
}