## sample from GEM process
rGEM <- 
  function(k = 500, alpha = 0.5, 
           plot = FALSE)
  {
    n <- k-1
    beta_raw <- rbeta(n, 1, alpha)
    beta_prod <- 
      beta_raw[-1] * 
      cumprod(1 - beta_raw)[-n]
    
    beta_out <- c(beta_raw[1], beta_prod)
    beta_out <- c(beta_out, 1-sum(beta_out))
    
    if (plot){
      x <- rnorm(k)
      plot(x, y = beta_out, type = "n")
      segments(x0 = x, y0 = 0, y1 = beta_out)
    }
    beta_out     
  }