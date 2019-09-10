# More accurately estimate average fragment length by summarizing peak overlapping on opposite strands
CHOPseq.PairedPeakDist<-function(query, subject, min.depth=0, max.depth=Inf, step=1, start=0, end=frag.len*2) {
#query, subject      GRanges objects, parameter passed to findOverlap() function, peaks on opposite strands
#min.depth   integer, minimum depth of query peaks
#max.depth   integer, maximum depth of query peaks
#step        integer, number of bases to move at each step
#start, end  integer, base range within which to search the maximum overlapping

library(chipseq);

# select peaks within given depth limits
depth<-elementMetadata(query)$depth;
query<-query[depth>=min.depth&depth<=max.depth];

names(query)<-1:length(query);
names(subject)<-1:length(subject);

# split by chromosomes
qr<-split(query, as.vector(seqnames(query)));
sb<-split(subject, as.vector(seqnames(subject)));
chr<-intersect(names(qr), names(sb));
qr<-qr[chr];
sb<-sb[chr];

# convert to IRanges objects
qr<-lapply(qr, ranges);
sb<-lapply(sb, ranges);

# set peak width to 1, located at the middle point
qr<-lapply(qr, function(qr) resize(qr, width=1, fix='center'));
sb<-lapply(sb, function(sb) resize(sb, width=1, fix='center'));

# set values to shift peaks
if (step==0) step<-1;
if (end>start&step<0) step<--1*step;
if (start>end&step>0) step<--1*step;
sh<-seq(start, end, step);

 getOlap<-function(qr, sb, gp, sh) {
 olap<-lapply(sh, function(sh, qr, sb, gp) findOverlaps(shift(qr, shift=sh), sb, maxgap=gp), qr=qr, sb=sb, gp=gp);
 olap<-lapply(olap, as.matrix); 
 distance<-rep(sh, sapply(olap, nrow));
 olap<-do.call(rbind, olap);
 olap[,1]<-as.numeric(names(qr)[olap[,1]]);
 olap[,2]<-as.numeric(names(sb)[olap[,2]]);
 cbind(distance, olap);
 }

##############################################################################################################################
olap<-lapply(chr, function(chr, qr, sb, sh, gp) getOlap(qr[[chr]], sb[[chr]], gp, sh), qr=qr, sb=sb, sh=sh, gp=abs(step)-1);
##############################################################################################################################

olap<-do.call(rbind, olap);
olap<-cbind(olap, elementMetadata(query)$depth[olap[,2]], elementMetadata(subject)$depth[olap[,3]]);
colnames(olap)<-c('Distance', 'Ind_query', 'Ind_subject', 'Depth_query', 'Depth_subject');

olap;
}
