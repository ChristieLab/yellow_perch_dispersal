######################################################
################# 05 October 2020 ####################
#################PCA's for GB and MB##################
######################################################

# Set the working directory 
setwd("~/Desktop/New_Final")

library(ggplot2)
library(adegenet)
library(hierfstat)
library(RColorBrewer)


#import genepop data and convert to genid object readable for adegenet
yp_indiv_final <- read.genepop("5781_noLD.gen", ncode = 3, quiet = TRUE)

all_names <- readLines("popvec_new.csv") #csv with pop assignments equal to n in pop. 
pop(yp_indiv_final) <- all_names #replaces incorrect pop names with new, correct names 
yp_all_final <- scaleGen(yp_indiv_final, NA.method = "mean") #turns our data into a matrix and replaces NA's

yp_indiv_final@pop 

dim(yp_all_final)
str(yp_all_final)
#run pca with prcomp
pca_all_final <- prcomp(yp_all_final, scale=FALSE) 
## plot pc1 and pc2
plot(pca_all_final$x[,1], pca_all_final$x[,2])

## make a scree plot
pca.var_final <- pca_all_final$sdev^2
pca.var.per_final <- round(pca.var_final/sum(pca.var_final)*100, 1)

#make pca results a data frame so we cna plot it with ggplot
pca.data_final <- data.frame(Sample=rownames(pca_all_final$x),
                       X=pca_all_final$x[,1],
                       Y=pca_all_final$x[,2])
pca.data_final

gbmb <- read.csv("GBMB.csv", header = FALSE) 

pca.data_final$basin = gbmb[,1]
pca.data_final

ggplot(data=pca.data_final, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = basin)) +
  xlab(paste("PC1 - ", pca.var.per_final[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per_final[2], "%", sep="")) +
  theme_bw() 

#okay, we want to do a pca with mb and gb samples separately as well. to do this we need to partition our pops
sep <- seppop(yp_indiv_final)
yp_mb<- repool(sep$CHI18, sep$WAK18, sep$NCH18, sep$SCH18, sep$GRH18, sep$GRH19, sep$MIC18, sep$MIC19, sep$CHE18, sep$CHX19, sep$MICYO, sep$MIL19, sep$STJ18, sep$SOH18, sep$SUT18, sep$NPT18, sep$NUB18, sep$MAN18, sep$LUD)
yp_gb <- repool(sep$BDN19, sep$BDNYO, sep$SGB18, sep$SGB19, sep$LBD19, sep$LBDYO, sep$MEN19)  


#let's start with green bay 
gb <- scaleGen(yp_gb, NA.method = "mean") #turns our data into a matrix and replaces NA's

#run pca with prcomp
pca_gb <- prcomp(gb, scale=FALSE) 
## plot pc1 and pc2
plot(pca_gb$x[,1], pca_gb$x[,2])

## make a scree plot
pca.var.gb <- pca_gb$sdev^2
pca.var.per.gb <- round(pca.var.gb/sum(pca.var.gb)*100, 1)

barplot(pca.var.per.gb, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#make pca results a data frame so we can plot it with ggplot
pca.data.gb <- data.frame(Sample=rownames(pca_gb$x),
                       X=pca_gb$x[,1],
                       Y=pca_gb$x[,2])

pca.data.gb
#let's add those populations names back in 
gb_pops <- read.csv("gb_pops.csv", header = TRUE) #read in pop name list for gb

pca.data.gb$pop = gb_pops$pop #add to existing dataframe
pca.data.gb$group = gb_pops$group

pca.data.gb #check

#plot by groups
ggplot(data=pca.data.gb, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = group), size = 3, alpha = .6 ) +
  scale_color_manual(breaks = c("BDN","LBD", "SGB"), values = c( "#f4a582","#d6604d","#67001f")) +
  xlab(paste("PC1 - ", pca.var.per.gb[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per.gb[2], "%", sep="")) +
  theme_bw() +
  theme(legend.title=element_blank())
  #theme(legend.position=c(.75, .75))   #if you wanna move your legend around


#plot by populations
ggplot(data=pca.data.gb, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = pop)) +
  scale_color_manual(breaks = c("BDN19", "BDNYO", "LBD19", "LBDYO", "MEN19", "SGB18", "SGB19"), values = c("palevioletred4", "palevioletred", "turquoise4", "paleturquoise", "yellow", "chartreuse4", "palegreen")) +
  xlab(paste("PC1 - ", pca.var.per.gb[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per.gb[2], "%", sep="")) +
  theme_bw() +
  theme(legend.title=element_blank())
#theme(legend.position=c(.75, .75))   #if you wanna move your legend around


#color numbers for chosen palette: 
brewer.pal(n = 10, "RdBu")



########################################################################### 
#now main basin
#let's try to remove the indivs skewing the pca

mb <- scaleGen(yp_mb, NA.method = "mean") #turns our data into a matrix and replaces NA's

remove <- c("CHI18_39", "IMC19_42", "MIL19_21", "MIL19_23", "GH19_6", "GH18_54", "MCY_15", "MCB18_25", "MCB18_21")
mb_fix <-mb[!row.names(mb) %in% remove, ]

#some of these alleles a fixed in the mb, giving 0 variance and throwing errors for prcomp. we want to remove them


#run pca with prcomp
pca_mb <- prcomp(mb_fix, scale=FALSE) 
## plot pc1 and pc2
plot(pca_mb$x[,1], pca_mb$x[,2])

## make a scree plot
pca.var.mb <- pca_mb$sdev^2
pca.var.per.mb <- round(pca.var.mb/sum(pca.var.mb)*100, 1)

barplot(pca.var.per.mb, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#make pca results a data frame so we can plot it with ggplot
pca.data.mb <- data.frame(Sample=rownames(pca_mb$x),
                          X=pca_mb$x[,1],
                          Y=pca_mb$x[,2])



#let's add those populations names back in 
mb_pops <- read.csv("mb_pops.csv", header = TRUE) #read in pop name list for mb

pca.data.mb$pop = mb_pops$pop #add to existing dataframe
pca.data.mb$region = mb_pops$region
pca.data.mb #check

#plot #########################################FIX
ggplot(data=pca.data.mb, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = region), size = 3, alpha = .6, ) +
  scale_color_manual(breaks = c("South", "North"), values = c("#2166AC", "#92C5DE")) +
  xlab(paste("PC1 - ", pca.var.per.mb[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per.mb[2], "%", sep="")) +
  theme_bw() +
  theme(legend.title=element_blank())

#theme(legend.position=c(.75, .75))   #if you wanna move your legend around

write.xlsx(pca.data.mb, file = "mb_pca_data_fix.xlsx", sheetName = "all", 
           col.names = TRUE, row.names = TRUE, append = FALSE)




