# Get segments of fixed length on a long vector (ex. depth on a chromosome) around a given set of indexes
# return a matrix with each row represents a segment
CHOPseq.GetSegments<-function(v, ind, before=0, after=before) {
# v; a numerical  vector or Rle object from which the segments are retrieved
# ind; indexes on the vector
# before; extended number of locations before the indexes
# after; extended number of locations after the indexes
library(chipseq);

v<-c(rep(NA, before), as.vector(v), rep(NA, after)); 
length<-before+after+1;


if (length(ind)>0) {
ss<-seqselect(v, start=ind, width=length); # return a vector with all the segments in concatenation
segments<-matrix(ss, nrow=length(ind), ncol=length, byrow=T);
}
else segments<-NA;

segments;
}

# Get segments of non-fixed length on a long vector (ex. depth on a chromosome) around a given set of indexes
# return a list of numeric vectors  with each element represents a segment
CHOPseq.GetSegments.2<-function(v, start, end, as.vector=F) {
# v; a numerical vector or Rle object from which the segments are retrieved
# start; Integer vector, start locations
# end; Integer vector, end locations, same length as 'start' 
# as.vector; if TRUE, convert Rle elements to vectors
library(chipseq);

start[start<1]<-1; start[start>length(v)]<-length(v); 
end[end<1]<-1; end[end>length(v)]<-length(v); 

x<-seqselect(v, start=start, end=end); # retrieve subset into a single vector
if (as.vector) x<-as.vector(x); 

size<-end-start+1; # size of ranges
factor<-rep(1:length(size), size); # splitting factor

segment<-split(x, factor);
segment;
}
