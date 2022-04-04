##################################################################
#########    Claire Schraidt 13 OCT 2020      #################### 
######### popstats and PW Fsts w/ bootstrap CI's #################
##################################################################

install.packages("adegenet", dep = TRUE)
install.packages("xlsx")
install.packages("hierfstat")
install.packages("sf") 

library(adegenet)
library(xlsx)
library(hierfstat)

# Set the working directory 
setwd("~/Desktop/New_Final")

#import genepop data and convert to genid object readable for adegenet
yp_indiv_hwe <- read.genepop("Hdplot_1snp_9302_927_final.gen", ncode = 3, quiet = TRUE)

#data should now appear as a GENIND object
yp_indiv_hwe

#look at pop names to make sure they've been converted correctly 
yp_indiv_hwe@pop

#they have not, so we need to rename them. 
all_names <- readLines("popvec_new.csv") #csv with pop assignments 
pop(yp_indiv_hwe) <- all_names #replaces incorrect pop names with new, correct names 


#now check the names
yp_indiv_hwe@pop  

#okay, PW fst. We need to convert out genind object to something usable by hierfstat
yp_hier_all <- genind2hierfstat(yp_indiv_hwe) #genind to hierfstat first...
str(yp_hier_all)


yp_pwfst <- pairwise.WCfst(yp_hier_all)
write.xlsx(LD_pwfst, file = "pw_fst_9302.xlsx", sheetName = "all", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

#use boot.ppfst to determine 95% CI's

pw_ci <- boot.ppfst(yp_hier_all, nboot=100,quant=c(0.025,0.975),diploid=TRUE)

write.xlsx(pw_ci[["ll"]], file = "ll.xlsx")
write.xlsx(pw_ci[["ul"]], file = "ul.xlsx")





#use hierfstat to caluclate popstats Ar, Hs, and Ho

ar_all <- allelic.richness(yp_hier_all,min.n=NULL,diploid=TRUE)
#export to excel and average pop values at each loci to get overall Ar
write.xlsx(ar_all$Ar, file = "Ar.xlsx")


stats<-basic.stats(yp_hier_all)

#export to excel and average pop values at each loci to get overall Ho and hs
write.xlsx(stats[["Ho"]], file = "Ho.xlsx")
write.xlsx(stats[["Hs"]], file = "Hs.xlsx")


