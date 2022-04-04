FST <- function(pop1, pop2) {



######################################################
## calculate FST between 2 populations (pop1 vs. pop2) - written by Xiaoshen Yin
######################################################

## read in genotype data for pop1 and pop2
#pop1<-read.table("C:/Users/Mark Christie/Downloads/pop1.txt",header=FALSE,sep="\t")
#pop2<-read.table("C:/Users/Mark Christie/Downloads/pop2.txt",header=FALSE,sep="\t")

## remoce the population column (i.e., column 1)
#pop1$V1<-NULL # xiaoshen used this when data was read in from file, not sure why
#pop2$V1<-NULL

#pop1 <- pop1[, -1]
#pop2 <- pop2[, -1]

## create dataframes to keep genotype data

popgeno1<-data.frame(matrix(,nrow=dim(pop1)[1],ncol=(dim(pop1)[2]/2)))
popgeno2<-data.frame(matrix(,nrow=dim(pop2)[1],ncol=(dim(pop2)[2]/2)))

for (i in 1:(dim(pop1)[2]/2)){
popgeno1[,i]<-paste0(pop1[,i*2-1],pop1[,i*2])
}

for (i in 1:(dim(pop2)[2]/2)){
popgeno2[,i]<-paste0(pop2[,i*2-1],pop2[,i*2])
}

## transpose genotyep dataframe to create dataframe with locus in row and sample in column
popgeno1_t<-as.data.frame(t(as.matrix(popgeno1)))
popgeno2_t<-as.data.frame(t(as.matrix(popgeno2)))

## calculate genotype and allele frequency
popgeno1_t$G11<-rowSums(popgeno1_t[,1:dim(pop1)[1]]=="11")
popgeno1_t$G12<-rowSums(popgeno1_t[,1:dim(pop1)[1]]=="12")
popgeno1_t$G21<-rowSums(popgeno1_t[,1:dim(pop1)[1]]=="21")
popgeno1_t$G22<-rowSums(popgeno1_t[,1:dim(pop1)[1]]=="22")
popgeno1_t$missing<-rowSums(popgeno1_t[,1:dim(pop1)[1]]=="-9-9")
popgeno1_t$gtotal<-rowSums(popgeno1_t[,c("G11","G12","G21","G22","missing")])
unique(popgeno1_t$gtotal)

popgeno2_t$G11<-rowSums(popgeno2_t[,1:dim(pop2)[1]]=="11")
popgeno2_t$G12<-rowSums(popgeno2_t[,1:dim(pop2)[1]]=="12")
popgeno2_t$G21<-rowSums(popgeno2_t[,1:dim(pop2)[1]]=="21")
popgeno2_t$G22<-rowSums(popgeno2_t[,1:dim(pop2)[1]]=="22")
popgeno2_t$missing<-rowSums(popgeno2_t[,1:dim(pop2)[1]]=="-9-9")
popgeno2_t$gtotal<-rowSums(popgeno2_t[,c("G11","G12","G21","G22","missing")])
unique(popgeno2_t$gtotal)

popgeno1_t$allele1<-(popgeno1_t$G11*2+popgeno1_t$G12+popgeno1_t$G21)/(2*popgeno1_t$G11+2*popgeno1_t$G12+2*popgeno1_t$G21+2*popgeno1_t$G22)
popgeno1_t$allele2<-(popgeno1_t$G22*2+popgeno1_t$G12+popgeno1_t$G21)/(2*popgeno1_t$G11+2*popgeno1_t$G12+2*popgeno1_t$G21+2*popgeno1_t$G22)
popgeno1_t$atotal<-popgeno1_t$allele1*(2*popgeno1_t$G11+2*popgeno1_t$G12+2*popgeno1_t$G21+2*popgeno1_t$G22)+popgeno1_t$allele2*(2*popgeno1_t$G11+2*popgeno1_t$G12+2*popgeno1_t$G21+2*popgeno1_t$G22)+popgeno1_t$missing*2
unique(popgeno1_t$atotal)

popgeno2_t$allele1<-(popgeno2_t$G11*2+popgeno2_t$G12+popgeno2_t$G21)/(2*popgeno2_t$G11+2*popgeno2_t$G12+2*popgeno2_t$G21+2*popgeno2_t$G22)
popgeno2_t$allele2<-(popgeno2_t$G22*2+popgeno2_t$G12+popgeno2_t$G21)/(2*popgeno2_t$G11+2*popgeno2_t$G12+2*popgeno2_t$G21+2*popgeno2_t$G22)
popgeno2_t$atotal<-popgeno2_t$allele1*(2*popgeno2_t$G11+2*popgeno2_t$G12+2*popgeno2_t$G21+2*popgeno2_t$G22)+popgeno2_t$allele2*(2*popgeno2_t$G11+2*popgeno2_t$G12+2*popgeno2_t$G21+2*popgeno2_t$G22)+popgeno2_t$missing*2
unique(popgeno2_t$atotal)

unique(popgeno1_t$allele1+popgeno1_t$allele2)
unique(popgeno2_t$allele1+popgeno2_t$allele2)

## create genotype and allele frequency dataframe

freq1<-data.frame(rownames(popgeno1_t),popgeno1_t$G11,popgeno1_t$G12,popgeno1_t$G21,popgeno1_t$G22,popgeno1_t$allele1,popgeno1_t$allele2)
colnames(freq1)<-c("locus","G11","G12","G21","G22","allele1","allele2")

freq2<-data.frame(rownames(popgeno2_t),popgeno2_t$G11,popgeno2_t$G12,popgeno2_t$G21,popgeno2_t$G22,popgeno2_t$allele1,popgeno2_t$allele2)
colnames(freq2)<-c("locus","G11","G12","G21","G22","allele1","allele2")

## merge two datasets to get shared loci for calculating FST
freq12<-merge(freq1,freq2,by.x="locus",by.y="locus")
colnames(freq12)<-c("locus","freq1.G11","freq1.G12","freq1.G21","freq1.G22","freq1.allele1","freq1.allele2","freq2.G11","freq2.G12","freq2.G21","freq2.G22","freq2.allele1","freq2.allele2")

# calculate FST
freq12$n1<-(freq12$freq1.G11+freq12$freq1.G12+freq12$freq1.G21+freq12$freq1.G22)
freq12$n2<-(freq12$freq2.G11+freq12$freq2.G12+freq12$freq2.G21+freq12$freq2.G22)

r<-2

freq12$nbar<-(freq12$n1+freq12$n2)/r
freq12$nc<-(r*freq12$nbar-((freq12$n1^2+freq12$n2^2)/(r*freq12$nbar)))/(r-1)

freq12$pbar1<-(freq12$n1*freq12$freq1.allele1+freq12$n2*freq12$freq2.allele1)/(r*freq12$nbar)
freq12$pbar2<-(freq12$n1*freq12$freq1.allele2+freq12$n2*freq12$freq2.allele2)/(r*freq12$nbar)

freq12$ssq1<-(freq12$n1*(freq12$freq1.allele1-freq12$pbar1)^2+freq12$n2*(freq12$freq2.allele1-freq12$pbar1)^2)/((r-1)*freq12$nbar)
freq12$ssq2<-(freq12$n1*(freq12$freq1.allele2-freq12$pbar2)^2+freq12$n2*(freq12$freq2.allele2-freq12$pbar2)^2)/((r-1)*freq12$nbar)

freq12$freq1.h1<-(freq12$freq1.G12+freq12$freq1.G21)/(freq12$freq1.G11+freq12$freq1.G12+freq12$freq1.G21+freq12$freq1.G22)
freq12$freq2.h1<-(freq12$freq2.G12+freq12$freq2.G21)/(freq12$freq2.G11+freq12$freq2.G12+freq12$freq2.G21+freq12$freq2.G22)

freq12$freq1.h2<-(freq12$freq1.G12+freq12$freq1.G21)/(freq12$freq1.G11+freq12$freq1.G12+freq12$freq1.G21+freq12$freq1.G22)
freq12$freq2.h2<-(freq12$freq2.G12+freq12$freq2.G21)/(freq12$freq2.G11+freq12$freq2.G12+freq12$freq2.G21+freq12$freq2.G22)

freq12$hbar1<-(freq12$n1*freq12$freq1.h1+freq12$n2*freq12$freq2.h1)/(r*freq12$nbar)
freq12$hbar2<-(freq12$n1*freq12$freq1.h2+freq12$n2*freq12$freq2.h2)/(r*freq12$nbar)

freq12$a1<-(freq12$nbar/freq12$nc)*(freq12$ssq1-(1/(freq12$nbar-1))*(freq12$pbar1*(1-freq12$pbar1)-(r-1)*freq12$ssq1/r-1/4*freq12$hbar1))
freq12$a2<-(freq12$nbar/freq12$nc)*(freq12$ssq2-(1/(freq12$nbar-1))*(freq12$pbar2*(1-freq12$pbar2)-(r-1)*freq12$ssq2/r-1/4*freq12$hbar2))

freq12$b1<-(freq12$nbar/(freq12$nbar-1))*(freq12$pbar1*(1-freq12$pbar1)-(r-1)*freq12$ssq1/r-(2*freq12$nbar-1)*freq12$hbar1/(4*freq12$nbar))
freq12$b2<-(freq12$nbar/(freq12$nbar-1))*(freq12$pbar2*(1-freq12$pbar2)-(r-1)*freq12$ssq2/r-(2*freq12$nbar-1)*freq12$hbar2/(4*freq12$nbar))

freq12$c1<-1/2*freq12$hbar1
freq12$c2<-1/2*freq12$hbar2

## consider a population that has a population structure of two levels
## one from the individual (I) to the subpopulation (S) and one from the subpopulation to the total (T)
## the total F, known here as FIT, can be partitioned into FIS (or f) and FST (or θ):

freq12$FST<-((freq12$a1+freq12$a2)/((freq12$a1+freq12$b1+freq12$c1)+(freq12$a2+freq12$b2+freq12$c2)))
freq12$FIT<-(((freq12$a1+freq12$b1)+(freq12$a2+freq12$b2))/((freq12$a1+freq12$b1+freq12$c1)+(freq12$a2+freq12$b2+freq12$c2)))
freq12$FIS<-((freq12$b1+freq12$b2)/((freq12$b1+freq12$c1)+(freq12$b2+freq12$c2)))

fsts <- freq12$FST
fst  <- mean(fsts, na.rm = TRUE)

#write.csv(freq12,file="FST_pop1v2.csv")
return(fst)

}

