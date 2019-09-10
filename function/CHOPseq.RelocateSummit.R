# Relocate summit position based on symmetry of peak shape
CHOPseq.RelocateSummit<-function(coverage, peaks, ws, min.ratio=2/3, max.ratio=3/2) {
# coverage  Rle, depth at all locations on a chromosome
# peaks     GRanges, peak information, the first 3 columns are start, end and depth
# ws        Integer, window size of the local region to search for a better summit
# min.ratio Numeric, cutoff of total depth at a position relative to the original summit position
# max.ratio Numeric, cutoff of total depth at a position relative to the original summit position
min.ratio<-min(min.ratio[1], 1);
max.ratio<-max(max.ratio[1], 1);

centers<-round(start(peaks)/2+end(peaks)/2); # initiate locations of peak summits

# retrieve depth information around peak centers
depth<-CHOPseq.GetSegments(coverage, centers, before=2*ws); 

# Total depth around each base
sums<-list(rowSums(depth[, 1:(1+2*ws)]))

for (i in 1:(2*ws)) sums[[i+1]]<-sums[[i]]-depth[,i]+depth[,i+1+2*ws];
sums<-do.call(cbind, sums);

    getDiff<-function(d1, d2, weight) {
    d<-abs(d1-d2);
    d<-t(t(d)*weight)/sum(weight);
    rowSums(d);
    }

diff<-lapply(1:(1+2*ws), function(i, depth, ws) getDiff(depth[, (i-1+ws):i], depth[, (i+1+ws):(i+2*ws)], sqrt(ws:1)), depth=depth, ws=ws);
diff<-do.call(cbind, diff);
diff<-diff/(sums+1);

ind<-lapply(1:nrow(diff), function(i, diff, sums, sum0) {
    x<-diff[i,];
    y<-sums[i,];
    x[y<(sum0[i]*min.ratio)|y>(sum0[i]*max.ratio)]<-max(x);
    
    # find local minimum
    a<-c(Inf, x[-length(x)]);
    b<-c(x[-1], Inf);
    ind<-which(x<a&x<b);

    if (length(ind)==0) ind<-1+ws;
    ind;
    }, diff=diff, sums=sums, sum0=sums[, 1+ws]);

# filter ind, so the depth at new location is at least 2/3 of the original
x<-rep(1:length(ind), sapply(ind, length));
y<-unlist(lapply(1:length(ind), function(i, ind, c, ws) ind[[i]]-ws-1+c[[i]], ind=ind, ws=ws, c=centers));
z<-rep(centers, sapply(ind, length));
r<-coverage[y]/coverage[z];
y[as.vector(r)<2/3]<-NA;
a<-split(y, x);
ind<-lapply(1:length(ind), function(i, ind, a) ind[[i]][!is.na(a[[i]])], ind=ind, a=a);

# select unique ind
ind<-sapply(ind, function(ind, ws) {
    # pick the nearest local minimum
    if (length(ind)>1) {
    d<-abs(ind-ws-1); # distance to center
    ind<-ind[d==min(d)];
    }    
    # If still 2, randomly pick one
    if (length(ind)==2) ind<-ind[sample(1:2, 1)];
    
    if (length(ind)==0) ind<-1+ws;
    ind;
    }, ws=ws);

elementMetadata(peaks)$adjusted.loc<-ind-ws-1+centers;
elementMetadata(peaks)$strand.diff<-sapply(1:length(ind), function(i, ind, diff) diff[i, ind[i]], ind=ind, diff=diff);
elementMetadata(peaks)$total.depth<-sapply(1:length(ind), function(i, ind, sums) sums[i, ind[i]], ind=ind, sums=sums);

peaks;
}