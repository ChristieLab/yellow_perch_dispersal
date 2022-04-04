#==================================================================================================#
# Script created by Mark Christie, contact at markchristie1500@gmail.com
# Script created in version R 4.0.2
# This script:  Generates model output for Larval Dispersal project
# Usage notes: Set all parameters below and then source this file
#==================================================================================================#
# Set working directory and output directory
# Directory where model.R and 'source' folder reside  

setwd("C:/Users/fishf/Dropbox/manuscripts/yellow perch rad-seq/model/larva3.1/")

base.directory <- getwd()
outdir <- paste(base.directory,"/output/",sep="")  # directory to save model output  
source(paste(base.directory, "/source/FunctionSourcer.R", sep = '')) #loads m.matrix, sources functions, and sets source directory

#==========================================================================================================#
n.replicates  <- 1    # number of replicates for each combination of parameters
replicates    <- Replicates(parameters, n.replicates)


for(i in 1:length(replicates[, 1])){


#samples <- 40
#for(s in 1:samples){
#  i = sample(1:length(replicates[, 1]), 1)

  n.generations      <<- as.numeric(as.character(replicates[i, 1]))
  adult.survival.var <<- as.numeric(as.character(replicates[i, 2]))
  myear              <<- as.character(replicates[i, 3])
  mmonth             <<- as.character(replicates[i, 4])
  mpld               <<- as.character(replicates[i, 5])
  
  model<- RunModel(n.generations, i=1)

  }


