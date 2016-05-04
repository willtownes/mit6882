library(matrixStats) #logSumExp functions

# generate forward probability 
forward_prob <- 
  function(x, z, lambda, pk, log_lik_func, msg, verbose){
    K <- nrow(pk)
    T <- nrow(x)
    
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
                 # handle z_{t-1} for t = 1
                 z_prev <- ifelse(t==1, z[1], z[t-1])
                 # compute
                 log(pk[z_prev, k]) +
                   log_lik_func(x[t, ], lambda[[k]]) + 
                   msg[k, t+1]
               }
        )
      prob_mat <- # normalize
        prob_mat - logSumExp(prob_mat)
      
      if (any(is.na(prob_mat)) | !all(is.finite(prob_mat))){
        stop("forward_prob: over/underflow prob value. Bad message values?")
      }
      prob[, t] <- exp(prob_mat)
    }
    
    # return
    prob
  }