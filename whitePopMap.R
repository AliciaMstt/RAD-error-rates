whitePopMap<-function(tsv, writedirectory, drop.rep, matinfo, PopOrder){
  # Function to create a list of the samples present in a desired matrix
  # and write a file that can be read by Stacks as a Population Map of desired samples to analyse with the populations program
  #
  #
  ## Variables:
  # tsv: the path to the file of the matrix of selected RADloci from the PostCleaning script
  # writedirectory: The directory where the PopMap file should be saved, whitout starting or ending with /
  # drop.rep: if TRUE replicates (sampled ending with _r or _ir) are discarded, if FALSE they are kept
  # matinfo: a cvs file with sample names in a column "sample", and population information in a column "Pop" 
  # PopOrder: numeric vector with desired order ti sort populations. Use to fit the PopKey desired order for your populations change line 43 (levels(x$Pop)….)
  
  # First load the matrix 
  final = read.delim(paste(tsv, sep = ""), header = T) 
  # Extact name of recovered samples
  samples <- colnames(final[, 4:ncol(final)])
  # Extract these samples from the matinfo
  x <- match(samples, matinfo$sample) #match labels with samples names
  x <- matinfo[x,] #create a new matrix with the samples from dat.d
  
  if(drop.rep == TRUE){
    # Drop replicates
    s <-grep("_r|_ir", x$sample, value= TRUE, invert= TRUE) # drop the samples ending with _r or _ir
    s <- match(s, x$sample) # select the rest
    x <- x[s,]
    # Keep only samples and population columns
    x <- x[, c(1,4)]
    # Change levels of Pop to integers
    levels(x$Pop) <- PopOrder
    # transform to vectors
    x$sample <- as.vector(x$sample)
    str(x)
    # Write it to a file with the name of the tvs matrix plus _PopMap_withrep (indicating reps are included)
    write.table(x, file= paste0(writedirectory, "/", basename(tsv), "_PopMap_norep.tsv"), sep = "\t",
      row.names =FALSE, col.names=FALSE,
      quote=FALSE)
  } else {
    # Don't drop replicates
    # Keep only samples and population columns
    x <- x[, c(1,4)]
    # Change levels of Pop to integers
    levels(x$Pop) <- PopOrder
    # transform to vectors
    x$sample <- as.vector(x$sample)
    str(x)
    # Write it to a file with the name of the tvs matrix plus _PopMap_norep (indicating reps are excluded)
    write.table(x, file= paste0(writedirectory, "/", basename(tsv), "_PopMap_withrep.tsv"), sep = "\t",
      row.names =FALSE, col.names=FALSE,
      quote=FALSE)   
  }
  
}
