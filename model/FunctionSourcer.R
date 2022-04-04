#FunctionSourcer <- function() {
# Set working directory, import packages, source functions, 
setwd(paste(base.directory,"/source/", sep = ''))    # set temp working directory 
#library(fields)  # used for calculation of pathogen dispersal distance matrix
#m.matrix <- read.table(paste(getwd(), "/connectivity.matrices/ocean_distance3.txt", sep = ''), header=FALSE)
#m.matrix <- read.table(paste(getwd(), "/connectivity.matrices/ocean_distance5.txt", sep = ''), header=FALSE)
#m.matrix <- read.table(paste(getwd(), "/connectivity.matrices/larval_connectivity_2014_1.txt", sep = ''), header=FALSE)


source(paste(getwd(), "/PopSetup.R", sep = ''))
source(paste(getwd(), "/Mortality.R", sep = ''))
source(paste(getwd(), "/NeededOffspring.R", sep = ''))
source(paste(getwd(), "/Pairs.R", sep = ''))
source(paste(getwd(), "/Repro.R", sep = ''))
source(paste(getwd(), "/Dispersal.R", sep = ''))
source(paste(getwd(), "/Parameters.R", sep = ''))
source(paste(getwd(), "/PlotIt.R", sep = ''))
source(paste(getwd(), "/Genotypes.R", sep = ''))
source(paste(getwd(), "/Output.R", sep = ''))
source(paste(getwd(), "/RunModel.R", sep = ''))
source(paste(getwd(), "/FST.R", sep = ''))
source(paste(getwd(), "/Replicates.R", sep = ''))
