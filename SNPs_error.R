SNP_error <- function(liSNPs, param){ # tsv is the SNP  matrix exported with export_sql.pl from stacks
  liSNPs = liSNPs
  
  ##### Function that 
  ## Uses a genlight object of SNPs
  ## to estimate SNP error rate
  ## replicate pairs should be labeled as _r or _ir at the end of the sample name

  ##### Variables:
  # liSNPs = genlight object (package adegenet) of SNPs (generated from e.g. a plink.raw to genlight object)
  # param = character string to label the stacks parameters or data used 
  
  
  require(adegenet)
### 1) Get recovered sample-replicate pairs
  # get samples names
  samples<-indNames(liSNPs)
  
  # get succesfull replicate-sample pairs
  reps <- grep("_r|_ir",samples,value=TRUE) # get the replicates (ending with _r or _ir)
  samps <- match(sub("_r|_ir","",reps),samples) # match against its sample (ie names w/o _r or _ir)
  samps<- samples[samps] # get names
  pairs<-cbind(reps,samps) # put them side by side 
  pairs<- pairs[rowSums(is.na(pairs)) < 1,] # remove rows that cointain NA to get only succesful replicate-sample pairs
  npairs <-nrow(data.frame(pairs))
  
### 2) Estimate SNP error rate 
y <- numeric(0)
for(i in 1:npairs) {
  
  srpair<-match(pairs[i,], samples) # 
  srpair <- liSNPs[srpair,] 
  indNames(srpair) # sample - replicate pair, first part of the name should be equal
  
  ### Estimate SNP-LOCI differences between each replicate pair  
  
  # estimate  number of loci and and missing loci (NAs)
  nloci <- nLoc(srpair) # number of total loci
  NAs<-NA.posi(srpair) # position of NAs in the sample-replicate pair
  
  
  # estimate number of NA in sample and replicate
  s<-NAs[1] # NAs in sample
  s<-s[[1]]
  
  r<- NAs[2] # NAs in replicate
  r<- r[[1]]
  
  # Find total missing loci between s and r togheter
  nMissLoc <- length(union(s,r))
  nMissLoc
  
  #### Estimate SNP differences between each replicate pair
  # handle the data as a matrix of SNPs
  mat <- as.matrix(liSNPs) 
  mat <- mat[indNames(srpair),]  # extract the replicate pairs
  
  # Estimate the number of loci excluding missing data, i.e. number of loci - the union of NA from both samples 
  n.loci.woNA <- nloci - nMissLoc   
  
  ## Estimate number of SNP mismatches 
  # The NA should be considered as missing data and not taken into the count.  
  
  # misSNP.cell = function to give a value of 1 if SNPs are different or 0 if not. For one locus, ignoring NA   
  misSNP.cell <- function(locus){ 
    all.r = mat[[1,locus]] # the SNP value for locus i in the replicate
    all.s = mat[[2,locus]] # the SNP value for locus i in the other replicate
    
    if(is.na(all.r) == FALSE & is.na(all.s) == FALSE) { # if both replicates do NOT have a missing locus 
      # prodceed to evaluate if they are the same or not
      SNP.diff <- all.s == all.r # evaluate if they are the same
      if(SNP.diff == TRUE) { return(0)  # give a value of 0 if they are the same
      } else { 
        misall = return(1)} # count as 1 difference if they are not the same
    } else {
      return(NA) 
    }
  }
  
  # Run in the whole mat
  misSNP.all <- mat.or.vec(1,ncol(mat))  # create a vector to store output
  for (i in 1:ncol(mat)) # run misSNP.cell in all the columns 
  {
    misSNP.all[i] <- misSNP.cell(i)
  }
  
  # and sum the result of each loci to get the total number of alelles that mismatch   
  SNP.mismatches<-sum(misSNP.all[1,], na.rm = TRUE)
  
  ## Estimate the SNP error rate as the number of SNP mismatches over the number of loci compared 
  SNP.error.rate <- SNP.mismatches/n.loci.woNA
  # And put results toghether by sample-replicate pair
  pair <- paste(indNames(srpair)[2],"-",indNames(srpair)[1], sep="") # to generate a name for the pair
  SNP.error.rate<-cbind(pair, SNP.error.rate)
  
  y<-rbind(y, SNP.error.rate)   
}

return(y)  
}