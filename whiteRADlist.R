whiteRADlist<-function(tsv, writedirectory){
  # Fuction to create a list of the CatalogID loci present in a desired matrix 
  # and write a file that can be read by Stacks as whitelist of desired loci
  # 
  ## Variables:
  # tsv the path to the file of the matrix of selected RADloci from the PostCleaning script
  # writedirectory The directory where the whitelist file should be saved, , whitout starting or ending /
  
  # First load the matrix 
  final = read.delim(paste(tsv, sep = ""), header = T) 
  # Extact name of loci (CatalogID)
  whitelist <- final$CatalogID
  # Write it to a file with the name of the tvs matrix plus _whitelist
  write(whitelist, file= paste0(writedirectory, "/", basename(tsv), "_whitelist.tsv"), ncolumns = 1)
}
