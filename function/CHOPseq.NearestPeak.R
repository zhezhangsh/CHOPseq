# Fine the nearest neighbor downstream of each peak. 
CHOPseq.NearestPeak<-function(fw, rv, tag.length, min.depth, ratio=.5) {
#fw,rv          data.frame or matrix, peaks found on forward and reverse streams of the same chromosome. The first 3 columns are start, end, and maximum depth of each peak
#tag.length		integer, average tag length
#min.depth      integer, minimum depth of peaks whose neighbor to be found
#ratio          numeric, the relative ratio of peak depth a neighboring peak must have. By default, if the peak has depth of 20, its neighbor must have depth of 10

len<-round(tag.length/2);
start<-c(fw[,1]-len, rv[,1]+len);
end<-c(fw[,2]-len, rv[,2]+len);
depth<-c(fw[,3], rv[,3]);
loc<-round(start/2+end/2);
strand<-rep(c(1, -1), c(nrow(fw), nrow(rv)));

peaks<-data.frame(adj.start=start, adj.end=end, adj.loc=loc, strand, depth);
peaks<-peaks[order(loc, strand), ];

	
	####################### get neighbor of peaks having given depth
	getNeighbors<-function(sortedPeaks, depth) {
	ind<-ind.nb<-which(sortedPeaks$depth==depth); 
	s<-sortedPeaks$strand[ind];
	
	ind.nb[s==1]<-ind[s==1]+1; # if peak on forward strand, find the neighbor downstream
	ind.nb[s==-1]<-ind[s==-1]-1; # if peak on reverse strand, find the neighbor upstream
	
	ind<-ind[ind.nb>0&ind.nb<=nrow(sortedPeaks)]; # not report neighbor if index out of bound
	ind.nb<-ind.nb[ind.nb>0&ind.nb<=nrow(sortedPeaks)];
	
	nb.dist<-abs(sortedPeaks$loc[ind]-sortedPeaks$loc[ind.nb]);
	same.strand<-sortedPeaks$strand[ind]*sortedPeaks$strand[ind.nb];
	
	cbind(sortedPeaks[ind, ], nb.depth=sortedPeaks$depth[ind.nb], nb.dist, same.strand);
	} #################### end of function
	
d<-as.numeric(names(table(depth[depth>=min.depth])));

nb<-lapply(d, function(d, peaks, ratio) getNeighbors(peaks[peaks$depth>=(d*ratio), ], d), peaks=peaks, ratio=ratio);

out<-do.call(rbind, nb);
}

