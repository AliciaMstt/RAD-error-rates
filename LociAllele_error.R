
LociAllele_error <- function(mat, param){ 
  ##### Function to estimate loci and allele error rates between replicate pairs based on tvs matrix (SNP.SNPs)
     ## replicate pairs should be labeled as _r or _ir at the end of the sample name
     ## outputs a data frame with the following columns for each replicate-pair found:
      ## "nloci" = total number of loci in the dataset 
      ## "nMissLoc" = number of missing loci for both samples of a replicate pair 
      ## "MissTotProp", = nMissLoc/nloci
      ## "shareMissLoc" = number of loci that are lost in both samples of  replicate pair
      ## "loci.mismatches" = number of loci that failed in the sample or replicate, but not in both
      ## "unshareMissLoc" = loci.mismatches/nMissLoc
      ## "loci.error.rate" =  loci.mismatches/nloci
      ## "n.loci.woNA" = number of loci that were NOT missing in both samples or a replicate pair
      ## "allele.mismatches" = number of alleles that mismatch between samples of a replicate pair 
      ## "allele.error.rate" allele.mismatches/n.loci.woNA
  
  ##### Variables:
  # mat =  dataframe containing a SNP matrix as the one exported with export_sql.pl from stacks
  # param = character string to label the stacks parameters or data used 
  
  final = mat
  
  ###### Function
  #### 1) convert data to a genlight object (adegenet package)
  ###  Prepare matrix to create genlight object
  ## function to transform the SNPs stored as text to an integer state
  snp2bin = function(loc){
    loc = as.character(unlist(loc))
    nsamp = length(loc)
    loc = as.factor(loc)
    nall = nlevels(loc)
    if(nall < 4){
      states = nchar(levels(loc))
      newlevs = c(0, 2, 1)
      newlevs = newlevs[order(states)]
      levels(loc) = newlevs
      as.numeric(as.character(loc))
    } else{
      rep(-9, nsamp)
    }
  }
  
  prepbin = function(i, ...){
    snp2bin(final[i, 4:ncol(final)]) #make sure loci data starts in col 3
  }
  
  ### convert to AA -> 0, aa -> 1, Aa -> 2 format
  recoded = sapply(1:nrow(final), prepbin)
  test = recoded[1, ] != -9
  test[is.na(test)] = T
  recoded = recoded[, test]
  colnames(recoded) = paste("locus", 1:ncol(recoded))
  rownames(recoded) = colnames(final[, 4:ncol(final)])
  
  recoded.l = NULL
  for(i in 1:nrow(recoded)){
    recoded.l[[i]] = recoded[ i,]
  }
  
  ### convert to genlight format
  require(adegenet)
  liSNPs <- new("genlight",  # stores data from many genotypes (samples) as binary SNPs
    gen= recoded.l, # the SNPs matrix coded as integers
    ploidy=2, 
    ind.names= colnames(final[, 4:ncol(final)])
  ) 
  
  #check data makes sense
  attributes(liSNPs)
  names(liSNPs)
  nInd(liSNPs) 
  nLoc(liSNPs)
  indNames(liSNPs)
  
  ### Get sample-replicate pairs
  # get samples names
  samples<-indNames(liSNPs)
  
  # get succesfull replicate-sample pairs
  reps <- grep("_r|_ir",samples,value=TRUE) # get the replicates (ending with _r or _ir)
  samps <- match(sub("_r|_ir","",reps),samples) # match against its sample (ie names w/o _r or _ir)
  samps<- samples[samps] # get names
  pairs<-cbind(reps,samps) # put them side by side 
  pairs<- pairs[rowSums(is.na(pairs)) < 1,] # remove rows that cointain NA to get only succesful replicate-sample pairs
  npairs <-nrow(data.frame(pairs))
  
  #### 2) Estimate loci, allele and error.rate  differences between each replicate pair
  
  # Create matrix to store data
  repliDiff<-matrix(ncol=12)
  colnames(repliDiff)<-c("parameter","pair", "nloci", "nMissLoc", "MissTotProp", "shareMissLoc", "loci.mismatches", "unshareMissLoc", "loci.error.rate", "n.loci.woNA", "allele.mismatches", "allele.error.rate")
  repliDiff<-repliDiff[-1,]
  
  for(i in 1:npairs) {
    
    srpair<-match(pairs[i,], samples) # 
    srpair <- liSNPs[srpair,] 
    indNames(srpair) # sample - replicate pair, first part of the name should be equal
    
    ### Estimate LOCI differences between each replicate pair  
    
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
    
    
    # Find proportion of missing loci over total number of loci
    MissTotProp  <- nMissLoc/nloci
    MissTotProp  
    
    # Find loci that are absent in both sample and replicate 
    #  (but that exist in other samples from the dataset)
    shareMissLoc <- length(intersect(s,r)) # elments common not sample and replicate
    shareMissLoc 
    
    # Find loci that failed in the sample or replicate, but not in both
    # this are a source of difference between replicates
    loci.mismatches <- union(setdiff(s,r), setdiff(r,s)) #elements in s not in r + elements in r not in s
    loci.mismatches <- length(loci.mismatches)
    loci.mismatches
    
    # Proportion of unshared loci out of the total missing loci
    unshareMissLoc <- loci.mismatches/nMissLoc
    unshareMissLoc  
    
    # Poroportion of unshared loci out of the total loci found: error rate per loci  
    loci.error.rate <- loci.mismatches/nloci
    loci.error.rate
    
    #### Estimate ALLELE differences between each replicate pair
    # handle the data as a matrix of alleles
    mat <- as.matrix(liSNPs) 
    mat <- mat[indNames(srpair),]  # extract the replicate pairs
    
    # Estimate the number of loci excluding missing data, i.e. number of loci - the union of NA from both samples 
    n.loci.woNA <- nloci - nMissLoc   
    
    ## Estimate number of allele mismatches 
    # The NA should be considered as missing data and not taken into the count.  
    
    # misalle.cell = function to give a value of 1 if alleles are different or 0 if not. For one locus, ignoring NA   
    misalle.cell <- function(locus){ 
      all.r = mat[[1,locus]] # the allele value for locus i in the replicate
      all.s = mat[[2,locus]] # the allele value for locus i in the other replicate
      
      if(is.na(all.r) == FALSE & is.na(all.s) == FALSE) { # if both replicates do NOT have a missing locus 
        # prodceed to evaluate if they are the same or not
        alle.diff <- all.s == all.r # evaluate if they are the same
        if(alle.diff == TRUE) { return(0)  # give a value of 0 if they are the same
        } else { 
          misall = return(1)} # count as 1 difference if they are not the same
      } else {
        return(NA) 
      }
    }
    
    # Run in the whole mat
    misalle.all <- mat.or.vec(1,ncol(mat))  # create a vector to store output
    for (i in 1:ncol(mat)) # run misalle.cell in all the columns 
    {
      misalle.all[i] <- misalle.cell(i)
    }
    # and sum the result of each loci to get the total number of alelles that mismatch   
    allele.mismatches<-sum(misalle.all[1,], na.rm = TRUE)
    
    ## Estimate the allele error rate as the number of allele mismatches over the number of loci compared 
    allele.error.rate <- allele.mismatches/n.loci.woNA
    
    ###### Put all B) results toghether by sample-replicate pair
    pair <- paste(indNames(srpair)[2],"-",indNames(srpair)[1], sep="") # to generate a name for the pair
    x<-cbind(param, pair, nloci, nMissLoc, MissTotProp, shareMissLoc, loci.mismatches, unshareMissLoc, loci.error.rate, n.loci.woNA, allele.mismatches, allele.error.rate)
    
    repliDiff<-rbind(repliDiff, x)
  }
     
return(as.data.frame(repliDiff, stringsAsFactors= FALSE))
}
