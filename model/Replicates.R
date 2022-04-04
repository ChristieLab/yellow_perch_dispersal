Replicates <- function(parameters, n.replicates) {
  # takes all *.set variables and creates the precise number of replicates for each variable
  
  #THINGS TO VARY: resistiance.cost, n.treated.streams, tfm strength, tfm timing - MAKE NOTE OF CURRENT PARAMETERS - seems to be working!
  
  
  n.generations       <- c(50, 100, 200)
  adult.survival.var  <- c(105)  # 105 is default; std.devation around mean survival each year
  
  # below is for selecting combinations of connectivity matrices
  #group  <- c("group2") # adding other groups needs debugging - currently 2014_6 (group1) == 2014_12 (group2) 
  years  <- c(2014, 2015, 2016, 2017, 2018, 2019, "2014-2019")
  months <- c(1, 2, 3, 4, 5, 6, "1-3", "4-6", "1-6", 7, 8, 9, 10, 11, 12, "7-9", "10-12", "7-12")
  pld    <- c(30, 40, 50)
  
  replicates <- expand.grid(n.generations, adult.survival.var, years, months, pld)
  replicates <- replicates[rep(seq_len(nrow(replicates)), n.replicates), ]
  
  return(replicates)
}  
