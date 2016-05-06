library(matrixStats) #logSumExp functions

# generate backward message
backward_message <- 
  function(y, x, z, theta, pk, 
           log_lik_func, verbose)
  {
    # OUTPUT:
    # msg: message for time 1 to T+1, d x (T + 1) dimension 
    
    T <- nrow(y)
    K <- nrow(pk)
    log_pk <- log(pk)
    
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
        log_lik_t <- 
          sapply(1:K, function(j)
            log_lik_func(x[t+1, ], x[t, ], theta[[j]])
          )
        msg[, t] <- colLogSumExps(log_pk + (log_lik_t + msg[, t+1]))
        #~~~~~~~~~~~~~~~~~~~~~
        # backup for old method
        # msg_mat <- # state-dish message matrix
        #   outer(
        #     1:K,1:K,
        #     Vectorize(
        #       function(j, k){
        #         log_pk[j, k] +
        #           log_lik_func(x[t+1, ], x[t, ], theta[[j]]) + 
        #           msg[j, t+1]
        #       }
        #     )
        #   ) 
        # 
        # msg[, t] <- # sum of state
        #   colLogSumExps(msg_mat)
      
      if (any(is.na(msg[, t])) | !all(is.finite(msg[, t]))){
        stop("backward_message: over/underflow message value. Bad lambda values?")
      }
    }
    # return
    msg
  }