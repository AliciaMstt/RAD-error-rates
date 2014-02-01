PairsDiff <- function(directory, tsv, plink, pat){ 
  ##### Function that 
  ## 1) converts data to a genlight object (adegenet package) to  
  ##  estimate loci, allele, and error rates between replicate pairs based on tvs matrix
  ## 2) Estimates SNP error rate based on a genlight object created from plink SNPs file
  ## 3) Visualize the Matrix of SNPs as a plot
  
  # directory: path to the directory where the SNP and the plinkf file are
  # tsv : is the SNP.SNPs  matrix exported with export_sql.pl from stacks and processed by the PostCleaning.r script
  # plink : is the plink file with SNP data corresponding to the tsv file
  # pat : is the pattern common to all the different files that would be compared. Eg: "JmonExplo", it 
    # would be used to identify the different parameters used (different values of m, M, n, etc)
  
  
  ## 
  directory = directory
  tsv = tsv
  plink = plink
  pat = pat
  
  
  # Get the name of the parameter used 
  pattern <- paste0("(.SNP.SNPs){0,1}(", pat, "_){0,1}(", pat, "){0,1}")
  param<-gsub(tsv, pattern=pattern, replacement="") 
  
  # get data
  final = read.delim(paste0(directory,tsv), header = T) 
  
  #### 1) Convert data to a genlight object (adegenet package) to estimate allele and loci error rates
  
  # Load custom function
  source(paste0(WD,"/bin/LociAllele_error.R"))
  
  # Run funtion to get allele and loci error rates
  repliDiff<-LociAllele_error(mat=final, param=param)  
  repliDiff
  
  # See results
  repliDiff
  
  ### 2) Estimates SNP error rate based on a genlight object created from plink SNPs file
  
  # Get the name of the parameter without deleting _
  pattern<- paste0("(.SNP.SNPs){0,1}(", pat, "){0,1}")
  param<-gsub(tsv, pattern=pattern, replacement="") 
  
  ## Create genlight object based on plink
  require(adegenet)
  liSNPs<- read.PLINK(file = paste0(directory,plink))
  
  # source SNP:error script
  source(paste0(WD,"/bin/SNPs_error.R"))
  
  # Transform param to id (do not include _)
  id<-gsub(tsv, pattern="(.SNP.SNPs){0,1}(JmonExplo_){0,1}(JmonExplo){0,1}", replacement="")
  
  # Run funtion to get SNP error rate
  y<-SNP_error(liSNPs=liSNPs, param=id)
  y<-as.data.frame(y)  
  
  # Add SNP.error.rate to the matrix of results
  repliDiff <-merge(repliDiff, y)
  
  # Save as object to compare against other params
  # Transform param to do not include _ again
  
  # Transform param to do not include _ 
  pattern<- paste0("(.SNP.SNPs){0,1}(", pat, "_){0,1}(", pat,"){0,1}")
  param<-gsub(tsv, pattern=pattern, replacement="")
  
  repliDiff.tsv <- paste("repliDiff.", param, sep="")
  assign(repliDiff.tsv, repliDiff, envir =.GlobalEnv)
  
  ##### 3) Visualize the Matrix as a plot
  # in this plot AA -> 0, aa -> 1, Aa -> 2 
  glPlot(liSNPs)
  title(paste("Distribution of SNPs by individual", "using", param))
  
  
}
