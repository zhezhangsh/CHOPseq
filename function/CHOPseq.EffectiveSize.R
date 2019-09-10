# get effective size of chromosomes based on their tag coverage
CHOPseq.EffectiveSize<-function(cov, gap.size=10000) {
# cov		list of Rle, each has the coverage information of a chormosome
# gap.size	integer, the minimum size of gaps void of tags, equal to the minimal size of assembly gaps in human genome
library(chipseq);
if (gap.size<=0) gap.size=10000
len<-sapply(cov, function(cov) {
	lengths<-cov@lengths;
	values<-cov@values;
	sum(lengths[lengths<10000|values>0]); }
);
names(len)<-names(cov);

len;
}

