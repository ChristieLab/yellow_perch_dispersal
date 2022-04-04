Mortality <- function(pops, death) {

  total.n <- length(pops[, 1])
  n.die   <- total.n * death  # number to remove
  #n.die   <- round(rnorm(1, mean = n.die, sd = adult.survival.var))  # variation in number that die
  pops    <- pops[-(sample(1:total.n, n.die, replace = FALSE)), ]
 
  return(pops)
}


