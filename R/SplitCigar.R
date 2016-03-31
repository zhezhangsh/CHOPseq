# Split cigar string
SplitCigar<-function(cigar, op=c('M', 'S', 'H', 'I', 'D')) {
  # cigar   Vector of cigar strings
  # op      Cigar operations
  
  require(GenomicAlignments);
  
  lst<-cigarToRleList(cigar); 
  
  val<-runValue(lst);
  len<-runLength(lst);
  ele<-elementLengths(len);
  
  val<-unlist(val, use.names=FALSE);
  len<-unlist(len, use.names=FALSE);
  ind<-rep(1:length(ele), ele); 
  
  v<-rep(0, length(cigar));
  n<-sapply(op, function(o) { print(val[1]); print(length(val)); 
    i<-which(val==o);
    s<-sapply(split(len[i], ind[i]), sum);
    v[as.integer(names(s))]<-s;
    v;
  });
  
  rownames(n)<-cigar;
  
  n; 
}