Repro = function(pops, allpairs, n.loci, id.un){
  # generation starting ID number
  OFFSPRING <- NULL
  for(n in 1:length(allpairs)) {
    
    pairs.n <- allpairs[[n]]
    pops.n  <- pops[pops[, 1] == n, ]
    #if(length(pairs.n) == 0) {pairs.n=1}  
    if(length(pairs.n) == 0) {next}  
    
    # randomly pair males and females
    size    <- min(table(pops.n[, 4]))
    males   <- sample(which(pops.n[, 4] == 1), size, replace = FALSE)
    females <- sample(which(pops.n[, 4] == 2), size, replace = FALSE)
    pairs   <- cbind(males, females)
    pairs   <- pairs[sample(1:length(pairs[, 1]), length(pairs.n), replace = TRUE), ] #randomly select needed number of pairs; replace = TRUE in case not enough pairs
    
    if(length(pairs) == 2) {pairs = cbind(pairs[1], pairs[2], pairs.n)} else {pairs   <- cbind(pairs, pairs.n)} # add number of offspring each pair will produce
    
    # set up all offspring matrix
    ALLOFFSPRING      <- matrix(nrow = sum(pairs.n), ncol = 7 + (n.loci * 2)) #IF VARIANCE EXCEEDS 0.25 EXTRA OFFSPRING, WILL GET ERROR
    ALLOFFSPRING[, 1] <- -9
    totaloffs         <- 0
  
    # generate offspring
    for(p in 1:nrow(pairs)) {
      f = pairs[p, 2]
      m = pairs[p, 1]
      noff = pairs[p, 3]
     
      
      if(noff > 0){
        #ids <- seq(from = gstartID + 1, to = gstartID + noff, by = 1)
        #population offpsring matrix
        offspring      = matrix(nrow=noff, ncol=5)
        offspring[,1]  = rep(n, noff)                                         # population residing
        offspring[,2]  = offspring[, 1]                                       # population of birth
        offspring[,3]  = rep(1, noff)                                         # Id
        offspring[,4]  = sample(c(1,2), noff, replace=TRUE, prob=c(0.5, 0.5)) # male=1, female=2
        offspring[,5]  = rep(0, noff)                                          # age
        
        # genotypes
        # prep parent genotypes
        fg = pops.n[f, -c(1:5)]
        mg = pops.n[m, -c(1:5)]
        fgid = pops.n[f, 3]
        mgid = pops.n[m, 3]
        
        # prep offspring genotype matrix
        offspringG = matrix(nrow=noff, ncol=(n.loci*2))
        
        # allele 1 positions
        positions = seq(1, (n.loci*2), 2)
        
        # randomly sample either position 1 or 2 (add 0 or 1) to starting position
        fallele  <- positions + sample(0:1, n.loci * noff, replace = TRUE)
        fallele2 <- fg[fallele]
        fallele3 <- matrix(fallele2, nrow = noff, ncol = n.loci, byrow = TRUE)
        
        mallele  <- positions + sample(0:1, n.loci * noff, replace = TRUE)
        mallele2 <- mg[mallele]
        mallele3 <- matrix(mallele2, nrow = noff, ncol = n.loci, byrow = TRUE)
        
        offspringG[, positions]     <- fallele3
        offspringG[, positions + 1] <- mallele3
        
        offspring    = cbind(offspring, offspringG)
        offspring    = cbind(offspring, fgid, mgid)
        #gstartID     = gstartID + nrow(offspring)
        
        m1 <- match(-9, ALLOFFSPRING[, 1])
        ALLOFFSPRING[m1:((m1+nrow(offspring))-1), ] <- offspring
      } 
    }
  
  OFFSPRING <- rbind(OFFSPRING, ALLOFFSPRING)
  }
  ids <- seq(from = id.un+1, to = (id.un+1) + (length(OFFSPRING[, 1])-1), by = 1) # add in unique ids
  OFFSPRING[, 3] <- ids
  return(OFFSPRING)
}

