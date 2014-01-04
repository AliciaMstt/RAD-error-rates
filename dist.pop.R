
#### Estimate distance between individuals of the same pop 
# storing it in a single data.frame

dists<-data.frame()
for (i in c("m3", "m4", "m10", "def")){

param <- i

### 1) Create genlight object 
  # Inport using the plink file exported from populations Stacks program
  # using whilelist file (i.e. loci in final matrix) and pop map excluding replicates
  require(adegenet)
  liSNPs<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_", param,"/Popsouts_Rselec/out.noreplicates/plink.raw"))
  
  # Change pop labels to Pop names
  levels(pop(liSNPs))<-c("Aj","An","Iz","Ma","Out", "Pe","Tl","To","Za")

### 2) Compute pairwise distances between individuals 
  
  ## separate the whole genlight object in smaller ones by creating blocks of loci. This facilitates computing
  blocks<- seploc(liSNPs, n.block=5) 
  class(blocks) #check if it is a list
  ## estimate distance matrix between individuals of each block
  D<- lapply(blocks, function(e) dist(as.matrix(e))) 
  names(D) #check names correspond to blocks
  # generate final general distance matrix by summing the distance matrixes
  Df<- Reduce("+", D) 
  
  ###  Compute distance matrix (euclidean)
  dat.d = dist(Df)
  
  ## add matinfo to dat.d matrix
  x <- match(labels(dat.d), matinfo$sample) #match labels with samples names
  x <- matinfo[x,] #create a new matrix with the samples from dat.d
  pop <- x$Pop #extract Pop info
  lane <-x$Lane #extract Lane info
    
### 3) Normalize distance matrix
  dat.d<- dat.d/max(dat.d)

### 4) Extract distances between indvs of each population
  dist.p<-data.frame()
  for (i in levels(pop(liSNPs))){
  # get the samples belonging to the i pop from the dad.d dist matrix
  s<-grep(i, labels(dat.d), value= TRUE, invert= FALSE) 
  # Transform dist matrix to matrix
  m<-as.matrix(dat.d)
  # select the samples to keep
  m<-m[s,s]
  # keep only lower triangle of the mat
  dist<-m[lower.tri(m)]
  # store data in the data.frame dists created externally 
  Pop <- rep(i,length(dist))  
  x<-cbind.data.frame(Pop, dist)
    
  dist.p<-rbind.data.frame(x,dist.p)
  }
  
  # add param info
  param<-rep(param,nrow(dist.p))
  x<-cbind.data.frame(param,dist.p)
  dists<-rbind.data.frame(x,dists)
}
