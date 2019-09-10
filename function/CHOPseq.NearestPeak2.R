# Fine the nearest neighbor downstream of each peak. 
CHOPseq.NearestPeak<-function(fw, rv, tag.ength, min.depth, ratio=.5) {
#fw,rv		data.frame or matrix, peaks found on forward and reverse streams of the same chromosome. The first 3 columns are start, end, and maximum depth of each peak
#tag.length		integer, average tag length
#min.depth	integer, minimum depth of peaks whose neighbor to be found
#ratio		numeric, the relative ratio of peak depth a neighboring peak must have. By default, if the peak has depth of 20, its neighbor must have depth of 10

start<-c(fw[,1], rv[,1]);
end<-c(fw[,2], rv[,2]);
depth<-c(fw[,3], rv[,3]);
loc<-round(start/2+end/2);
strand<-rep(c(1, -1), c(nrow(fw), nrow(rv)));

peaks<-data.frame(start, end, loc, strand, depth);
peaks<-peaks[order(loc, strand), ];

	
	####################### get neighbor of peaks having given depth
	getNeighbors<-function(sortedPeaks, depth) {
	ind<-which(sortedPeaks$depth==depth); 
	ind.bf<-ind-1; if (ind[1]==1) ind.bf=ind[1]+1;
	ind.af<-ind+1;	if (ind[length(ind)]==nrow(sortedPeaks)) ind.af[length(ind)]<-ind.bf[length(ind)]; 
	
	ind.nb<-ind.bf
	ind.nb[abs(sortedPeaks$loc[ind.af])<abs(sortedPeaks$loc[ind.bf])]<-ind.af[abs(sortedPeaks$loc[ind.af])<abs(sortedPeaks$loc[ind.bf])];

	nb.dist<-sortedPeaks$strand[ind]*(sortedPeaks$loc[ind]-sortedPeaks$loc[ind.nb]);
	same.strand<-sortedPeaks$strand[ind]*sortedPeaks$strand[ind.nb];
	
	cbind(sortedPeaks[ind, ], nb.depth=sortedPeaks$depth[ind.nb], nb.dist, same.strand);
	} #################### end of function
	
d<-as.numeric(names(table(depth[depth>=min.depth])));

nb<-lapply(d, function(d, peaks, ratio) getNeighbors(peaks[peaks$depth>=(d*ratio), ], d), peaks=peaks, ratio=ratio);

do.call(rbind, nb);
}

