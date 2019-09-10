# get the average of peak density in different chromosomal regions
CHOPseq.PeakDensity<-function(peaks, regions) {
# peaks     GRanges, locations of peaks
# regions   GRanges, different types of chromosomal regions

library(IRanges);

# The total size of each type of regions, calculate from scratch if not exists
total.size<-metadata(regions)$total.size
if (is.null(total.size)) {
region.type<-elementMetadata(regions)$Region
types<-unique(region.type);
total.size<-sapply(types, function(x, y, z) {a<-z[y==x]; a<-coverage(a); a[a>1]<-1; sum(sum(a));}, y=region.type, z=regions)
}

# names of region types
types<-names(total.size);

strand(peaks)<-'*';
strand(regions)<-'*';

olap<-as.matrix(findOverlaps(peaks, regions));
x<-split(olap[,1], elementMetadata(regions)$Region[olap[,2]]);
counts<-sapply(x, function(x) length(unique(x)));

types<-intersect(names(counts), names(total.size));
counts[types]/total.size[types];
}