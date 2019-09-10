# Pair peaks on opposite strands
CHOPseq.PairPeaks<-function(peak.fw, peak.rv, maxgap, min.depth=0) {
# peak.fw, peak.rv    GRanges objects, Peaks identified from the 2 strands, required peak unique ID and metadata fields: 'depth', 'adjusted.loc'
# maxgap            Parameter passed to findOverlaps() function
# min.depth         Integer, minimum depth of the query peaks
fw<-GRanges(seqnames=as.vector(seqnames(peak.fw)), ranges=IRanges(start=elementMetadata(peak.fw)[, 'adjusted.loc'], end=elementMetadata(peak.fw)[, 'adjusted.loc']));
rv<-GRanges(seqnames=as.vector(seqnames(peak.rv)), ranges=IRanges(start=elementMetadata(peak.rv)[, 'adjusted.loc'], end=elementMetadata(peak.rv)[, 'adjusted.loc']));

over.fw<-findOverlaps(fw, rv, maxgap=maxgap);
over.fw<-as.matrix(over.fw);

over.rv<-findOverlaps(rv, fw, maxgap=maxgap);
over.rv<-as.matrix(over.rv);

    getResults<-function(over, peak1, peak2, min.depth) {
    chr<-as.vector(seqnames(peak1))[over[,1]];
    d1<-elementMetadata(peak1)$depth[over[,1]];
    d2<-elementMetadata(peak2)$depth[over[,2]];
    l1<-elementMetadata(peak1)$adjusted.loc[over[,1]];
    l2<-elementMetadata(peak2)$adjusted.loc[over[,2]];
    x<-data.frame(Chr=chr, Query.loc=l1, Query.depth=d1, Subject.loc=l2, Subject.depth=d2);
    x[x[,3]>=min.depth,];
    }

map.fw<-getResults(over.fw, peak.fw, peak.rv, min.depth);
map.rv<-getResults(over.rv, peak.rv, peak.fw, min.depth);
map.rv[, c(2, 4)]<--1*map.rv[, c(2,4)];
map<-rbind(map.fw, map.rv);

id<-paste(as.vector(map[,1]), map[,2], map[,4], sep='_');
split(1:nrow(map), id)->x;
x<-sapply(x, function(x) x[1]);
map<-map[sort(x),];
ratio<-map[,5]/map[,3];
map<-map[ratio>=.2&ratio<=4,];
rownames(map)<-1:nrow(map);

d<-map[,4]-map[,2];
x<-density(d);
adj0<-round(x$x[which(x$y==max(x$y))]);
map<-data.frame(map, Distance.adj=d-adj0);

# map<-split(map[, -1], as.vector(map[,1]));

    adj<-function(d) {
    x<-density(d);
    round(x$x[which(x$y==max(x$y))]);
    }

map
}
