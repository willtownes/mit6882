library(matrixStats) #logSumExp functions

# generate backward message
backward_message <- 
  function(x, lambda, pk, log_lik_func, verbose){
    K <- nrow(pk)
    T <- nrow(x)
    
    # initialize message, a state*time matrix size K*T
    msg <- matrix(NaN, nrow = K, ncol = T+1)
    msg[, T+1] <- 0
    
    # compute backward message (log scale)
    if (verbose){
      cat("\n Backward Message Passing:\n")
      pb <- txtProgressBar(0, T, style = 3)
    }
    for (t in T:1){
      if (verbose) setTxtProgressBar(pb, t)
      msg_mat <- # state-dish message matrix
        outer(
          1:K,1:K,
          Vectorize(
            function(j, k){
              log(pk[j, k]) +
                log_lik_func(x[t, ], lambda[[j]]) + 
                msg[j, t+1]
            }
          )
        )
      msg[, t] <- # sum of state
        colLogSumExps(msg_mat)
      
      if (any(is.na(msg[, t])) | !all(is.finite(msg[, t]))){
        stop("backward_message: over/underflow message value. Bad lambda values?")
      }
    }
    
    # return
    msg
  }