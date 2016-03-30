# Split cigar string
SplitCigar<-function(cigar, op=c('M', 'S', 'H', 'I', 'D')) {
  # cigar   Vector of cigar strings
  # op      Cigar operations
  
  require(GenomicAlignments);
  
  lst<-cigarToRleList(cigar); 
  
  val<-lapply(lst, runValue);
  len<-lapply(lst, runLength); 
  ele<-elementLengths(len);
  
  val<-unlist(val, use.names=FALSE);
  len<-unlist(len, use.naems=FALSE);
  ind<-unlist(1:length(cnt), ele); 
  
  v<-rep(0, length(cigar));
  n<-sapply(op, function(o) {
    i<-val==o;
    s<-sapply(split(len[i], ind[i]), sum);
    v[as.integer(names(s))]<-s;
    v;
  });
  
  n; 
}