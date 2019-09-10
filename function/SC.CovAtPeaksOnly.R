# Keep regions (islands) overlapping to peaks; set the background regions as 0
SC.CovAtPeaksOnly<-function(peaks, coverage, min.depth=1) {
# peaks         GRanges, locations of peaks
# coverage	SimpleRleList, read depth info of all genome
# min.depth     Integer, minimum depth to show in bedGraph

islands<-slice(coverage, lower=min.depth);
islands<-GRanges(seqnames=space(islands), ranges=IRanges(unlist(start(islands)), unlist(end(islands))));
seqlengths(islands)<-sapply(coverage, length);

n<-countOverlaps(islands, peaks);
islands<-islands[n>0];
c<-coverage(islands); # 0 if overlap to peaks, 1 otherwise
    
c*coverage; # regions out of peaks will be 0; regions in peaks will be unchanged
}   
