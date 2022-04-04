PopSetup <- function(n.pops, k.adults) {

  total.n <- n.pops * k.adults

  pops <- sort(rep(1:n.pops, k.adults))
  #ids  <- rep(1:k.adults, n.pops)
  ids  <- pops
  pid  <- 1:total.n
  sex  <- sample(1:2, total.n, replace = TRUE)  # create sexes in 50:50 sex ratio
  year <- 0

  population <- cbind(pops, ids, pid, sex, year)
  
  return(population)
}


