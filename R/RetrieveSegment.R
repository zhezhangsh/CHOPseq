# retrieve depth information of given segments with fixed width from all chromosomes
RetrieveSegmentList<-function(values, ranges, before=0, after=before) {
  # values      A list of numeric vectors or Rles
  # ranges      an GRanges objects defining the ranges of data to be selected
  if (is.null(names(ranges))) names(ranges)<-1:length(ranges);
  
  inds<-start(resize(ranges, 1));
  names(inds)<-1:length(inds);
  inds<-split(inds, as.vector(seqnames(ranges)));
  inds<-inds[names(inds) %in% names(values)];
  
  ids<-unlist(lapply(inds, names), use.names=FALSE);
  
  v<-lapply(names(inds), function(chr) RetrieveSegment(values[[chr]], inds[[chr]], before=before, after=after));
  v<-do.call('rbind', v);
  rownames(v)<-ids;
  
  out<-matrix(nr=length(ranges), nc=ncol(v));
  rownames(out)<-names(ranges);
  out[as.integer(rownames(v)), ]<-v;
  
  out;
}

# Get segments of fixed length on a long vector (ex. depth on a chromosome) around a given set of indexes
# return a matrix with each row represents a segment
RetrieveSegment<-function(v, ind, before=0, after=before) {
  # v; a numerical  vector or Rle object from which the segments are retrieved
  # ind; indexes on the vector
  # before; extended number of locations before the indexes
  # after; extended number of locations after the indexes
  library(chipseq);
  
  v<-c(Rle(rep(NA, before)), v, Rle(rep(NA, after))); 
  length<-before+after+1;
  
  if (length(ind)>0) {
    #ss<-seqselect(v, start=ind, width=length); # return a vector with all the segments in concatenation
    ss<-v[IRanges(start=ind, width=length)]; 
    segments<-matrix(as.vector(ss), nrow=length(ind), ncol=length, byrow=T);
  }
  else segments<-NA;
  
  segments;
}
