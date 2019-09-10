########################
# scan the base coverage on a chromsome and identify peaks with a minimum depth
CHOPseq.PeakCaller.2<-function(coverage, min.depth, ws, strand="*", merge=T, verbose=F) {
# coverage  list of Rle, each element gives the depth info. of a chromosome
# min.depth Integer, minimum depth at peak summit
# ws        Integer, window size of 2 nearest peaks
# strand    Character, indicate the strand
# merge	    Boolean, whether to merge nearby peaks with the same height as the last step
library(chipseq);

l<-length(coverage);
	
index<-lapply(coverage, function(coverage, depth) which(coverage>=depth), depth=min.depth);

# initiate peaks
peaks<-lapply(1:l, function(i, c, index) CHOPseq.InitiatePeaks(c[[i]], index[[i]]), c=coverage, index=index);
if (verbose) print(paste("Peaks initiated, N =", sum(sapply(peaks, nrow))));

# remove small peaks that have higher depth within 50 bases around
peaks<-lapply(1:l, function(i, c, p) CHOPseq.RemoveSmallPeaks(c[[i]], p[[i]]), c=coverage, p=peaks);
if (verbose) print(paste("Small peaks removed, N =", sum(sapply(peaks, nrow))));

# create GRanges
chr<-Rle(names(coverage), sapply(peaks, nrow)); 
peaks<-do.call(rbind, peaks);
if (identical(strand, '+') | identical(strand, 1)) str<-rep('+', nrow(peaks))
else if (identical(strand, '-') | identical(strand, -1)) str<-rep('-', nrow(peaks))
else str<-rep('*', nrow(peaks));
gr<-GRanges(seqnames=chr, ranges=IRanges(start=peaks[,1], end=peaks[,2]), strand=str);

seqlengths(gr)<-sapply(coverage, length);
elementMetadata(gr)<-peaks[-(1:2)]; 

# Adjust peak location
peaks<-split(gr, seqnames(gr));
peaks<-peaks[intersect(names(peaks), names(coverage))];
coverage<-coverage[names(peaks)];
peaks<-lapply(1:length(peaks), function(i, c, p, ws) CHOPseq.RelocateSummit(c[[i]], p[[i]], ws=ws), c=coverage, p=peaks, ws=ws);
if (verbose) print('Peak summits re-located');

    # Whether peaks are local best
    isBest<-function(peaks, ws) {
    meta<-elementMetadata(peaks); 
    a<-IRanges(start=meta[,'adjusted.loc'], end=meta[,'adjusted.loc']);
    
    over<-findOverlaps(a, a, maxgap=ws);
    over<-as.matrix(over); 
    #over<-over[over[,1]!=over[,2],];
    
    x<-meta[over[,1],];
    y<-meta[over[,2],];
    z<-rep(1, nrow(over));
    z[y[,'depth']>x[,'depth']]<-0;
    z[y[,'depth']==x[,'depth']&y[,'strand.diff']<x[,'strand.diff']]<-0;
    z[y[,'depth']==x[,'depth']&y[,'strand.diff']==x[,'strand.diff']&y[,'total.depth']>x[,'total.depth']]<-0;
    
    isBest<-split(z, over[,1]);
    b<-sapply(isBest, min);
    c<-rep(T, length(b));
    c[b==0]<-F;
    c;
    }
    
peaks<-lapply(peaks, function(peaks, ws) peaks[isBest(peaks, ws)], ws=ws);
print(paste('Peaks filtered, N =', sum(sapply(peaks, length))));

gr<-do.call('c', peaks);
names(gr)<-1:length(gr);
gr;
}
