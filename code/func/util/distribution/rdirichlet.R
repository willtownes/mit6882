# sample from dirichlet distribution
# credit: http://goo.gl/j7BGf7
rdirichlet <- function(a) {
  y <- rgamma(length(a), a, 1)
  return(y / sum(y))
}