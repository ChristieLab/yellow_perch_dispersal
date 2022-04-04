Output <- function(pops, n) {
  

  #calculate  Weir and Cockerham Fst with base R
  n.cols <- ncol(pops)
  gpops  <- pops[, -c(2:5,n.cols-1, n.cols )]
  pops1  <- gpops[, 1]
  gtypes <- gpops[, -1]
  
  head(gpops)
  popids <- unique(pops[, 1])
  
  
  pairs <- t(combn(popids, 2))
  OUT  <- NULL
  for(n in 1:nrow(pairs)){
    pair <- pairs[n, ]
    pop1 <- gpops[which(gpops[, 1] == pair[1]), ]
    pop2 <- gpops[which(gpops[, 1] == pair[2]), ]
    pop1 <- pop1[, -1]
    pop2 <- pop2[, -1]
    
    fsts <- FST(pop1, pop2)
    out  <- cbind(t(pair), fsts)
    OUT  <- rbind(OUT, out)
  }  
    
  model <- OUT
  

  # compare model fst to empirical fst: empirical_wcfsts.txt was calculated with xiaoshens script (as are model fsts)
  empirical <- read.table("empirical_wcfsts.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
  head(empirical)
  emp_ids   <- paste(empirical[, 1], empirical[, 2], sep="_")
  empirical <- cbind(emp_ids, empirical)

  head(model)
  model_ids1  <- paste(model[, 1], model[, 2], sep="_")
  model_ids2  <- paste(model[, 2], model[, 1], sep="_")
  model      <- cbind(model_ids1, model_ids2, model)

  # begin analyses
  head(empirical)
  head(model)
  
  # find similar (need to merge because 1_3 can equal 3_1)
  m1 <- match(empirical[, 1], model[, 1]) # should be NAs
  m2 <- match(empirical[, 1], model[, 2])
  m3 <- which(is.na(m1) == TRUE)
  
  m1[m3] <- m2[m3]
  
  model_reduced <- model[m1, ]
  
  out <- cbind(empirical, model_reduced) # first two columns are real data (with wrong row/col numbers), but correct)
  m5  <- which(is.na(out[, 11]) == TRUE) # remove NAs
  out <- out[-m5, ]
   #out <- out[, -8] # remove extra name column
  
  #cols: empirical afd =6, model afd = 11, empirical fst = 15
  #plot(out[, 6], out[, 11]) # empirical afd vs model afd
  
  # output1 = all points========================================================#

   # below is incorrect, want to predict empirical fst (code left as reminder)
  #dat <- lm(as.numeric(out[, 11]) ~ as.numeric(out[, 6]))
  #summary(dat)
  
  dat <- lm(as.numeric(out[, 6]) ~ as.numeric(out[, 11]))
  #summary(dat)
  #plot(out[, 11], out[, 6])
  #abline(dat)
  
  
  # out1 = all points, empirical FST vs model FST
  output1 <- cbind(t(summary(dat)$coefficients[2, c(1, 4)]), summary(dat)$r.squared)
  cor1    <- suppressWarnings(cor.test(as.numeric(out[, 6]), as.numeric(out[, 11]), method = "spearman"))
  sh1     <- shapiro.test(as.numeric(out[, 6]))
  sh2     <- shapiro.test(as.numeric(out[, 11]))
  output1 <- cbind(output1,  cor1[[4]], sh1[[1]], sh2[[1]]) 
  
  # output3  = main basin ========================================================#
  
  # remove greenbay # matches to grids 1,2,3 on grid regions map
  m3 <- which(out[, 9] == 3)
  m2 <- which(out[, 9] == 2)
  m1 <- which(out[, 9] == 1)
  m32 <- which(out[, 10] == 3)
  m22 <- which(out[, 10] == 2)
  m12 <- which(out[, 10] == 1)
  
  m4   <- c(m3, m2, m1, m32, m22, m12)
  out2 <- out[-m4, ]

  
  dat     <- lm(as.numeric(out2[, 6]) ~ as.numeric(out2[, 11]))
  #summary(dat)
  #plot(out2[, 11], out2[, 6])
  #abline(dat)
  
  output3 <- cbind(t(summary(dat)$coefficients[2, c(1, 4)]), summary(dat)$r.squared)
  cor1    <- suppressWarnings(cor.test(as.numeric(out2[, 6]), as.numeric(out2[, 11]), method = "spearman"))
  sh1     <- shapiro.test(as.numeric(out2[, 6]))
  sh2     <- shapiro.test(as.numeric(out2[, 11]))
  output3 <- cbind(output3,  cor1[[4]], sh1[[1]], sh2[[1]]) 
  
  
  # output5  = main basin removing GRH and LUD========================================================#

  m1   <- which(out2[, 5] == "GRH19")
  m2   <- which(out2[, 4] == "GRH19")
  m3   <- c(m1, m2)
  out2 <- out2[-m3, ]
  
  m1   <- which(out2[, 5] == "LUD18")
  m2   <- which(out2[, 4] == "LUD18")
  m3   <- c(m1, m2)
  out2 <- out2[-m3, ]
  
  dat <- lm(as.numeric(out2[, 6]) ~ as.numeric(out2[, 11]))
  #summary(dat)
  #plot(out2[, 11], out2[, 6])
  #abline(dat)
  
  output5 <- cbind(t(summary(dat)$coefficients[2, c(1, 4)]), summary(dat)$r.squared)
  cor1    <- suppressWarnings(cor.test(as.numeric(out2[, 6]), as.numeric(out2[, 11]), method = "spearman"))
  sh1     <- shapiro.test(as.numeric(out2[, 6]))
  sh2     <- shapiro.test(as.numeric(out2[, 11]))
  dmean   <- mean(as.numeric(dist(out2[, 11])))
  dmax    <- max(as.numeric(dist(out2[, 11])))
  output5 <- cbind(output5,  cor1[[4]], sh1[[1]], sh2[[1]], dmean, dmax)
  
  
  # output5  = main basin removing GRH and LUD========================================================#
  
  
  m4 <- which(out[, 2] == 1 | out[, 3] == 1)
  m5 <- which(out[, 2] == 2 | out[, 3] == 2)
  m6 <- which(out[, 2] == 3 | out[, 3] == 3)
  m6 <- c(m4, m5, m6)
  
  out3 <- out[m6, ]
  
  
  m3 <- which(out3[, 2] == 3 & out3[, 3] == 3)
  m2 <- which(out3[, 2] == 2 & out3[, 3] == 3)
  m1 <- which(out3[, 2] == 1 & out3[, 3] == 3)
  
  m4 <- which(out3[, 2] == 3 & out3[, 3] == 2)
  m5 <- which(out3[, 2] == 2 & out3[, 3] == 2)
  m6 <- which(out3[, 2] == 1 & out3[, 3] == 2)
  
  m7 <- which(out3[, 2] == 3 & out3[, 3] == 1)
  m8 <- which(out3[, 2] == 2 & out3[, 3] == 1)
  m9 <- which(out3[, 2] == 1 & out3[, 3] == 1)
  
  m10 <- c(m1, m2, m3, m4, m5, m6, m7, m8, m9)
  out4 <- out3[-m10, ]
  #plot(out4[, 11], out4[, 6])
  
  
  dat <- lm(as.numeric(out4[, 6]) ~ as.numeric(out4[, 11]))
  #summary(dat)
  #plot(out4[, 11], out4[, 6])
  #abline(dat)
  
  output7 <- cbind(t(summary(dat)$coefficients[2, c(1, 4)]), summary(dat)$r.squared)
  cor1    <- suppressWarnings(cor.test(as.numeric(out4[, 6]), as.numeric(out4[, 11]), method = "spearman"))
  sh1     <- shapiro.test(as.numeric(out4[, 6]))
  sh2     <- shapiro.test(as.numeric(out4[, 11]))
  output7 <- cbind(output7, cor1[[4]], sh1[[1]], sh2[[1]])
  
  
## prepeare data for export ===================================================================================#
  
  output <- cbind(output1, output3, output5, output7)
  
  
  output <- cbind(i,  n.generations, adult.survival.var, myear, mmonth, mpld, output)
  write.table(output,  paste(outdir,"output.txt", sep=""), row.names = FALSE, col.names = FALSE, sep="\t", append = TRUE)
  # output 1 = all populations
  # output 3 = all main basin
  # output 5 = all main basin with luddington and grandhaven19 dropped
  # output 7 = main basin vs greenbay only
  # col headers = estimate of slope, slope p-val, correlation, spearmans correlation, shapiro empirical, shapiro model
  # output5 has 2 additional columns (max gap space and average gap space), to further examine bimodal data in model data only
  }
