# Calculate the correlation between coverage between forward and reverse strands of each chromosome and summarize all chromosomes
SC.StrandCorr<-function(gr, shift=0, resize=1, round=FALSE, summary.method='weighted.mean') {
# gr			GRanges, tags
# resize  If greater than 0, resize reads to given number
# round   If TRUE, round depth higher than 1 to 1
# shift		Integer vector, if 0, no shifting 
# summary.method	Character, the method to summarize results from multiple chromosomes
#					"weighted.mean", default, mean of chromosomes weighted by number of reads
#					"mean", simple mean, not recommended
#					"z", Fisher's r to z conversion. See www.statisticshell.com/meta.pdf or Field, A.P. (2001) Psychological Methods, 6(2), 161-180 for details

library(GenomicRanges); 

if (resize[1]>0) gr<-resize(gr, resize[1]);
counts<-table(seqnames(gr)); # number of reads on each chromosome

cov.fw<-coverage(gr[strand(gr)=='+']); # forward strands
cov.rv<-coverage(gr[strand(gr)=='-']); # reverse strands
cov.all<-coverage(gr); 
chr.names<-names(cov.all);

# remove chromosome(s) having no coverage on neither strand
max<-pmax(sapply(cov.fw, max), sapply(cov.rv, max));
cov.fw<-cov.fw[max>0];
cov.rv<-cov.rv[max>0];
counts<-counts[names(cov.fw)];

if (round) {
# to reduce the weight of ranges with many tags (repetitive regions in many cases)
cov.fw<-lapply(cov.fw, function(c) {c@values[c@values>1]<-1; c;});
cov.rv<-lapply(cov.rv, function(c) {c@values[c@values>1]<-1; c;});
}
############################################## calculate Pearson's correlation coefficients
corr<-lapply(chr.names, function(chr, fw, rv, shift) shiftApply(shift, fw[[chr]], rv[[chr]], FUN=cor), shift=shift, fw=cov.fw, rv=cov.rv); 
##############################################

names(corr)<-chr.names;
corr<-as.data.frame(corr); 
rownames(corr)<-shift;

# the effective length is the total length of chromosome without 10+kb regions void of tags, usually correspond to assembly gaps
#used to calculated weighted correlation
#if (setequal(names(effective.size), chr.names)) effective.len<-effective.size[chr.names]
#else {
#effective.len<-sapply(cov.all, function(c) {l<-c@lengths; v<-c@values; sum(l[l<=10000|v>0]);} );
#effective.len<-effective.len[chr.names];
#}

All<-apply(corr, 1, FUN=function(corr, n, method) SC.MetaCorr(corr, n, method), n=counts, method=summary.method);

corr<-cbind(All, corr); 

corr;
} 


