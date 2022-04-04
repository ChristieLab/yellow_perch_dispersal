RunModel <- function(n.generations, i) {
  i <<- i

  # hosts columns are as follows: individual, species, genotype, stage(tadpole=1;adult=2), infected(no=1,yes=2), xposition, yposition

  pops <- PopSetup(n.pops, k.adults) # Cols: population, id in population, unique id, sex, age
  pops <- Genotypes(pops, n.loci, n.alleles) # creates geneotypes
  pops <- cbind(pops, 0, 0) # create progenitor mother and father ids
  
  
  for(n in 1:n.generations){  
    
    # selecting appropriate connectivity matrix year
    if(nchar(myear) > 4) {
    years <- seq(as.numeric(substr(myear, 1, 4)), as.numeric(substr(myear, 6, 9))) 
    tyear <- as.character(sample(years, 1))} else {tyear <- myear} 
    
    # selecting appropriate connectivity matrix month
    if(nchar(mmonth) > 2) {
      month  <- seq(as.numeric(substr(mmonth, 1, 1)), as.numeric(substr(mmonth, 3, 4)))
      
      if(nchar(mmonth) == 5) {
        month  <- seq(as.numeric(substr(mmonth, 1, 2)), as.numeric(substr(mmonth, 4, 5)))}
    
      tmonth <- as.character(sample(month, 1))} else {tmonth <- mmonth}  

    
    dispersal.matrix <- read.table(  paste(getwd(), "/matrices/", "/", tyear, "_", tmonth, "_", mpld, "/", "larval_connectivity.txt", sep = '') , header=FALSE, skip = 2)
    
    id.un    <- max(pops[, 3])  # max unique id
    pops     <- Mortality(pops, death=death.sim)
    noffs    <- NeededOffspring(pops, k.adults=k.adults, dispersal.matrix, adult.survival.var=adult.survival.var)
    allpairs <- Pairs(noffs, g.shape, g.rate)
    repro    <- Repro(pops, allpairs, n.loci, id.un)
    pops[, 5]<- pops[, 5] + 1 # age remaining adults before offspring are introduced
    pops     <- Dispersal(pops, repro, noffs)
    #print(cbind(myear, mmonth, tyear, tmonth)) 
    
   #}
    
    if(n == n.generations){
      out      <- Output(pops, n)
    }
  }
  
}



