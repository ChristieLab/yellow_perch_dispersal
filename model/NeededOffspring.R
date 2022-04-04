NeededOffspring <- function(pops, k.adults, dispersal.matrix, adult.survival.var) {
  
  pop.n <- data.frame(table(pops[, 1])) # calculate current pop size
  #pop.n <- pop.n[order(as.numeric(as.character(pop.n[, 1]))), ] # orders pops by name (very important for Dispersal.R)
  needed.n <- as.numeric(k.adults - pop.n[, 2]) # calculated number of needed offspsring in each population
  
  OUT <- NULL
  for(n in 1:length(needed.n)){
    out <- round(rnorm(1, mean = needed.n[n], sd = adult.survival.var))  # variation in number that die 
    OUT <- c(OUT, out)
  }
  needed.n <- OUT # above introduces density independence
  needed.n[needed.n <= 0] <- 10
  #lapply(needed.n, rnorm, n=1, mean=needed.n, sd = adult.survival.var)
  
  zeros <- which(colSums(dispersal.matrix) == 0)
  
  # if statement is for matrices with zero columns
  if(length(zeros) > 0) {
    for(l in 1:length(zeros)){
        dispersal.matrix[zeros[l], zeros[l]] <- 0.00001}
        }  # add very small number to zeros columns
  
  # function to use multinomial to use dispersal matrix to calculate n.offspring from each population
  offsample <- function(x) {rmultinom(1:n.pops, needed.n[x], prob = dispersal.matrix[, x])}
  offs <- sapply(1:length(needed.n), offsample)
  
  # if statement is for matrices with zero columns
  if(length(zeros) > 0) {
    for(l in 1:length(zeros)){
      offs[zeros[l], zeros[l]] <- 10} # just add 1 offspring to column with all zeros
  }  # add very small number to zeros columns
  
  
  
  return(offs)

  }


