Genotypes <- function(pops, n.loci, n.alleles) {
  n.gtypes <- length(pops[, 3])
  
  a = c(1:n.alleles)
  b = rev(a)
  c = b/a
  d = sum(c)
  freqs = c/d # for snps will be 0.8,0.2
  
  lowestallele = 1
  alleles3 = cbind(seq(lowestallele, lowestallele+n.alleles-1, 1), freqs)
  #alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   # table allele frequencies
  homos=(alleles3[,2])^2                                                          # create homozygote allele frequencies
  homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
  hets=t(combn(alleles3[,2],2))                                                   # create heterozygote allele frequencies
  hetfreq=2*(hets[,1]*hets[,2])
  hetvals=t(combn(as.character(alleles3[,1]),2))                                  # create heterozygote allele names
  hets2=cbind(hetvals,hetfreq)
  gfreqs=rbind(hets2,homos2)                                                      # combine hets and homos and create genotypes
  
  n=n.gtypes * 10                                                                 # sample size of all simulated genotypes (* 10 for a buffer)
  gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            # create genotypes(by coloumn, 1 for each allele)
  gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
  gtypes=cbind(gfreqs1,gfreqs2)
  gtypes=gtypes[sample(1:length(gtypes[,1]),replace=FALSE),]                          # shuffle genotypes, consider omitting?
  gtypes <- cbind(as.numeric(gtypes[, 1]), as.numeric(gtypes[, 2]))
  
  OUT <- NULL  
  for (l in 1:n.loci){ 
    sg1=gtypes[sample(1:length(gtypes[,1]), n, replace=TRUE),]  #replace changed from FALSE to TRUE HERE ON 3.38.16 #WHAT DOES IT DO?
    OUT <- cbind(OUT, sg1)
  }
  #OUT <- cbind(1:length(OUT[, 1]), OUT)
  output    <- OUT[sample(1:length(OUT[, 1]), n.gtypes, replace = FALSE), ]
  output    <- cbind(1:n.gtypes, output) # i is the current iteration
  pops      <- cbind(pops, output)
  pops      <- pops[, -6]
  return(pops)
}
