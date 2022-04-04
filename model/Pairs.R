Pairs = function(noffs, g.shape, g.rate){
  
  #hist(rgamma(n.offs*10, shape = 1, rate = 0.5), breaks = 50) # histogram of distribution
  noffs2 = rowSums(noffs)
  
  all.pairs <- list()
  for(n in 1:length(noffs2)) {
  
    n.offs <- noffs2[n]
    if(n.offs == 0) {next} # if 0 offspring, go to next
    pair.rs <- rgamma(n.offs*10, shape = g.shape, rate = g.rate)  # sample from distribution
    pair.rs <- round(pair.rs) # round to get whole numbers
    pair.rs    
    # this whole next part is to get us exactly n.offs offspring
    # basic procdure is to get slightly less than n.offs (ocasionally will get exactly n.offs) and randomly add pairs with 1 offspring until reach n.offs
    pair.rs <- c(0, pair.rs) # fixes a bug where if no 0s the next line returns an empty vector
    pair.rs <- pair.rs[-which(pair.rs == 0)] # remove all 0s
    pair.rs2 <- cumsum(pair.rs)  # take cumlative sum of pair.rs
    
    #if first value is larger than the number of offspring needed in total, have one pair produce all offspring
    if(pair.rs2[1]>n.offs){
      all.pairs[[n]] = n.offs
    } else {
      
    closest  <- which(abs(pair.rs2-n.offs) == min(abs(pair.rs2-n.offs))) # the first x individuals in pair.rs are cloest to n.offs
    if(length(closest) > 1) {closest <- closest[1]}
    
    n.sampled <- pair.rs2[closest] # are we > or < n.offs
    if(n.sampled > n.offs) {closest <- closest-1} # if > n.offs take one less so that we are always less than n.offs
    
    rs <- pair.rs[1:closest] # take all pairs
    
    add.off <- n.offs-sum(rs)  # how many extra offspring do we need to add
    pairs   <- c(rs, rep(1, add.off))
    all.pairs[[n]] <- pairs
  }
}

  return(all.pairs)
  #sum(rs) # shoould equal n.offs
  #hist(rs, breaks = 20)
  
}






