## sample from GEM process
rMulti <- 
  function(n, beta){
    rmultinom(n, 1, beta) %>% 
      apply(2, function(x) which(x == 1))
  }