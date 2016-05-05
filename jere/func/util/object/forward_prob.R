library(matrixStats) #logSumExp functions

# generate forward probability 
forward_prob <- 
  function(y, x, z, msg,
           theta, pk, 
           log_lik_func, verbose)
  {
    T <- nrow(y)
    K <- nrow(pk)
    lpk <- log(pk)
    
    # initialize message, a state*time matrix size K*T
    prob <- matrix(NaN, nrow = K, ncol = T)
    
    # compute forward probability 
    if (verbose){
      cat("\n Forward Message Passing:\n")
      pb <- txtProgressBar(0, T, style = 3)
    }
    
    for (t in 1:T){
      if (verbose) setTxtProgressBar(pb, t)
      prob_mat <- # state-dish message matrix
        sapply(1:K,
               function(k){
                 # assume z[0] = z[1]
                 z_prev <- ifelse(t==1, z[1], z[t-1])
                 # compute (NOTE: x[t, ] refer to x at t-1)
                 log(pk[z_prev, k]) +
                   log_lik_func(x[t+1, ], x[t, ], theta[[k]]) + 
                   msg[k, t+1]
               }
        )
      
      prob_mat <- # normalize
        prob_mat - logSumExp(prob_mat)
      
      if (any(is.na(prob_mat)) | !all(is.finite(max(prob_mat, 0)))){
        stop("forward_prob: over/underflow prob value. Bad message values?")
      }
      prob[, t] <- exp(prob_mat)
    }
    
    # return
    prob
  }