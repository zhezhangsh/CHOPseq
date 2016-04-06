# Split cigar string
SplitCigar<-function(cigar, op=c('M', 'S', 'H', 'I', 'D')) {
  # cigar   Vector of cigar strings
  # op      Cigar operations
  
  require(GenomicAlignments);
  require(S4Vectors);
  
  lst<-cigarToRleList(cigar); 
  
  val<-runValue(lst);
  len<-runLength(lst);
  ele<-elementLengths(len);
  
  val<-BiocGenerics::unlist(val, use.names=FALSE);
  len<-BiocGenerics::unlist(len, use.names=FALSE);
  ind<-rep(1:length(ele), ele); 
  
  v<-rep(0, length(cigar));
  n<-sapply(op, function(o) { 
    i<-which(val==o);
    if (length(i)>0) {
      s<-sapply(split(len[i], ind[i]), sum);
      v[as.integer(names(s))]<-s;
    }
    v
  });
  if (class(n)!='matrix') n<-matrix(n, nc=length(op)); 
  
  rownames(n)<-cigar;
  colnames(n)<-op;
  
  n; 
}

# Split a very long vector of CIGAR string
SplitLongCigar<-function(cigar, op=c('M', 'S', 'H', 'I', 'D'), len=10^6, cores=2) {
  # cigar     A very long vector of CIGAR strings
  # len       The number of CIGAR strings to process per batch
  # cores     Number of cores for parallel processing
  
  require(GenomicAlignments);
  require(S4Vectors);
  require(parallel);
  require(CHOPseq); 
  
  if (length(cigar) <= len) SplitCigar(cigar, op) else {
    ind<-rep(1:ceiling(length(cigar)/len), each=len)[1:length(cigar)];
    cigar<-split(cigar, ind); 
    n<-parallel::mclapply(cigar, SplitCigar, mc.cores=max(2, cores)); 
    do.call('rbind', n); 
  }
  
}