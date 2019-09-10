# Summarize read statistics at each peak
CHOPseq.PeakStats<-function(peaks, gr, extend.before=0, extend.after=extend.before, by.strand=T) {
# peaks: a list of data frame, each corresponds to a chromosome; the first 2 fields of data frame are the start and end of peaks
# gr: GRanges object, total reads of a sample
# extend.before: extend the peak boundary by given number of bases before peaks, the suggested value is the extended read length -1
# extend.after: extend the peaks boundary by given number of bases after peaks
# by.strand: whether get the summarization of strands separately;
# return a list of data frame each corresponds to a chromosome if by.strand is FALSE; 
# return 3 such lists if by.strand is TRUE for both strands and their summation
library(chipseq);

if (is.null(names(peaks))) names(peaks)<-1:length(peaks); 

if (by.strand) {
gr.fw<-gr[strand(gr)=='+']; 
gr.rv<-gr[strand(gr)=='-'];

fw<-CHOPseq.TagCountAtPeak(peaks, gr.fw, extend.before, 0);
rv<-CHOPseq.TagCountAtPeak(peaks, gr.rv, 0, extend.after);
#x<-intersect(names(fw), names(rv));
fw<-fw[names(peaks)];
rv<-rv[names(peaks)];

fw.total<-elementMetadata(fw)$Tag.Total;
fw.unique<-elementMetadata(fw)$Tag.Unique;
fw.adj<-sqrt(fw.total)*sqrt(fw.unique);

rv.total<-elementMetadata(rv)$Tag.Total;
rv.unique<-elementMetadata(rv)$Tag.Unique;
rv.adj<-sqrt(rv.total)*sqrt(rv.unique);

adj<-2*sqrt(fw.adj)*sqrt(rv.adj);
percent.fw<-fw.adj/(fw.adj+rv.adj);

elementMetadata(peaks)$Adj.Tag.Count<-adj;
elementMetadata(peaks)$Percent.Forward<-percent.fw;
elementMetadata(peaks)$Total.Tag.Forward<-fw.total;
elementMetadata(peaks)$Unique.Tag.Forward<-fw.unique;
elementMetadata(peaks)$Total.Tag.Reverse<-rv.total;
elementMetadata(peaks)$Unique.Tag.Reverse<-rv.unique;

peaks;
}

else CHOPseq.TagCountAtPeak(peaks, gr, extend.before, extend.after);
}


# Calculated average depth of peaks between 2 locations
CHOPseq.DepthAtPeaks<-function(coverage, peaks) {
# coverage	Rle list, with base coverage of whole chromosomes
# peaks		GRanges, peaks

chr<-names(coverage); 
peaks<-split(peaks, seqnames(peaks));

chr<-intersect(names(peaks), chr);
if (length(chr)==0) stop("No common chromosomes");

cov<-coverage[chr];
peaks<-lapply(peaks[chr], ranges);

depth<-lapply(chr, FUN=function(chr, cov, peaks) aggregate(cov[[chr]], start=start(peaks[[chr]]), end=end(peaks[[chr]]), FUN=mean), cov=cov, peaks=peaks);

names(depth)<-chr;
depth;
}

# Count number of reads whose start positions are within given ranges
CHOPseq.TagCountAtPeak<-function(peaks, gr, extend.before=0, extend.after=extend.before) {
# peaks		GRanges object, peak locations
# gr		GRanges object, location of all tags 
# extend	Integer, number of bases to extend the boundaries of the peaks
# return the same peaks GRanges object with 3 extra metadata elements:
	# total number of reads in each region, 
	# number of bases having reads, 
	# and the maximum number of reads at a single base

cov<-coverage(resize(gr, width=1)); # make depth at each base equal to number of reads starting at the base

if (is.null(names(peaks))) names(peaks)<-1:length(peaks);
peaks<-split(peaks, seqnames(peaks));

chr<-intersect(names(cov), names(peaks));
if (length(chr)==0) stop("No common chromosomes");

cov<-cov[chr];
peaks<-peaks[chr];

start<-lapply(peaks, function(x) {y<-start(x)-extend.before; names(y)<-names(x); y;} );
end<-lapply(peaks, function(x) end(x)+extend.after);

 getCounts<-function(cov, start, end) { # coverage, start, and end on a chromosome
 if (length(start)>0) {
 counts<-CHOPseq.GetSegments.2(cov, start, end, as.vector=T); # when the returned segments are vectors, following 'lapply' is much faster
 counts<-lapply(counts, function(x) x[x>0]); 
 names(counts)<-names(start); 
 counts
 }
 }; # end of getCounts

counts<-lapply(chr, function(chr, cov, start, end) getCounts(cov[[chr]], start[[chr]], end[[chr]]), cov=cov, start=start, end=end); 

peaks<-unlist(peaks, use.names=F);
counts<-unlist(counts, recursive=F);
counts<-counts[names(peaks)];

elementMetadata(peaks)$Tag.Total<-sapply(counts, sum);
elementMetadata(peaks)$Tag.Unique<-sapply(counts, length);

peaks;
}



