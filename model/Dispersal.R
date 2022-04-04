Dispersal <- function(pops, repro, noffs) {

  LARVA <- NULL
  for (i in 1:length(unique(repro[, 1]))){
    larva <- repro[repro[, 1] == i, ] #isolate all offspring produced in population 1
    if(length(larva) == (ncol(pops))) {larva = t(larva)} # fix for if only 1 larval
    dest  <- rep(1:ncol(noffs), noffs[i, ]) #get new population ids for all offspring produced in population i going to population x
    larva[, 1] <- dest # change pop id of all individuals
    LARVA <- rbind(LARVA, larva) # rbind to new object; can't overwrite because old i's become new i's and will keep getting overwritten
  }
  
  pops <- rbind(pops, LARVA) # add offspring to pops object
  pops <- pops[order(as.numeric(as.character(pops[, 1]))), ] # orders pops by name to keep indiviudals together
  return(pops)
}


