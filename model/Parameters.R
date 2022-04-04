#Note: this is the default file; all variables are manimpulated in Replicates.R
# If a parameter has a value of NA; it is varied directly in Replicates.R

#Input parameters
parameters <- list()          # all parameters must be added to this list (see code below section break : ======)

#Timing parameters
#Population parameters
n.pops  	<- 40          # number of poulations, must be mirrored in connectivity matrices
death.sim <- 0.2         # What percentage of adults die each year 
g.shape 	<- 0.5         # Controls the shape of gamma distribution used for variance in reproductive success
g.rate  	<- 0.1         # Control the shape of gamma distribution
#hist(rgamma(1000, shape = g.shape, rate = g.rate), breaks = 50, xlab="Number of Offspring per Pair") 

#within pop.parameters
k.adults            <- 1000 # what is carrying capcity of each local population
adult.survival.var  <- 105  # std.devation around mean survival each year

n.generations  <- 50  # number of years to run full model (large pops and reduced mortality)
start.geneflow <- 1  # year to start gene flow (don't want it for testing)

n.loci    <- 100
n.alleles <- 2 #optimized for 2, likely wont work without testing for n > 2


# Add all parameters to list===============================================================#
parameters[["n.pops"]]    <- n.pops
parameters[["death.sim"]] <- death.sim
parameters[["g.shape"]]   <- g.shape
parameters[["g.rate"]]    <- g.rate          

parameters[["k.adults "]]          <- k.adults 
parameters[["adult.survival.var"]] <- adult.survival.var
parameters[["n.generations"]]      <- n.generations
parameters[["start.geneflow"]]     <- start.geneflow

parameters[["n.loci "]]   <- n.loci 
parameters[["n.alleles"]] <- n.alleles


