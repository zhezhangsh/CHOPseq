# Calculate the correlation between coverage between forward and reverse strands separately for tags with different maping quality scores
# Optionally, run a stepwise process which adds scores one by one and tracks the resulting strand correlation
CHOPseq.MapqCorr<-function(gr, scores=-1, chr=NA, stepwise=T, stepwise.order='highest', stepwise.next.best=F, verbose=F) {
# gr					GRanges, tags, must have a metadata field called 'mapq', accessed by elementMetadata(gr)$mapq
# scores				integer vector, the mapq scores to be included in the analysis
# chr                   character vector, names of chromosomes to use
# stepwise				boolean, if TRUE, also run a stepwise process to add mapq scores one by one and get the resulting strand correlation
# stepwise.order		character, if stepwise if TRUE, how to order the mapq scores. *if stepwise.next.best is TRUE, only the 1st score is selected based on this order*
#						'highest', default, score with the highest values
#						'best', score has the best strand correlation 
#						'largest', score has the largest number of tags
# stepwise.next.best	boolean, if TURE, the stepwise process will add the mapq score improve the strand correlation the most at each step
# WARNING!!! THIS IS A VERY TIME-CONSUMING OPTION, USE WITH CAUTION. 
library(chipseq);

if (!('mapq' %in% names(elementMetadata(gr)))) stop("\'mapq\' field missed from the tag metadata"); 
if (identical(chr, NA)) chr=names(seqlengths(gr));

gr<-gr[as.vector(seqnames(gr)) %in% chr]; 
if (length(gr)==0) stop("No tags from given chromosome(s)");

mapq<-elementMetadata(gr)$mapq;
scores.gr<-as.numeric(names(table(mapq))); 
# limit tags to existing and specified mapq scores 
if (!identical(scores, -1)) scores<-intersect(scores.gr, as.numeric(scores)) else scores<-scores.gr;
if (length(scores)<1) stop('No valid mapq scores'); 

by.mapq<-lapply(scores, function(score, gr, mapq) gr[mapq==score], gr=gr, mapq=mapq);
names(by.mapq)<-scores;
score.count<-sapply(by.mapq, length); 
if (verbose) print('Finished splitting by mapq');

    ############### round higher coverage to 1, for correlation calculation, it reduces the weight of tag-rich regions
	round1<-function(cov) { # a list of Rle with chromosome coverage
    lapply(cov, function(cov) {values<-cov@values; values[values>1]<-1; values->cov@values; cov;});
    }##############
    
cov<-lapply(by.mapq, function(tags) list(forward=round1(coverage(tags[strand(tags)=='+'])[chr]), reverse=round1(coverage(tags[strand(tags)=='-'])[chr]))); 
rm(by.mapq); # release memory
if (verbose) print('Finished converting to depth');
	
# get effective size of each chromosome, 10kb+ regions void of tags are not included
effective.len<-CHOPseq.EffectiveSize(coverage(gr)[chr]); 
if (verbose) print('Finished calculating effective chromosome size'); 

	#########################################
	getCorr<-function(cov, len) { # cov, 2 element list for 2 strands
		# No depth more than 1. Do this to avoid regions with higher depth (many in repetitive regions) are given too much weight in correlation
		# cov<-lapply(cov, function(cov) lapply(cov, function(cov) {values<-cov@values; values[values>1]<-1; values->cov@values; cov;}));
		
		r<-sapply(1:length(cov[[1]]), function(i, x, y) cor(x[[i]], y[[i]]), x=cov[[1]], y=cov[[2]]); # a numeric vector with correlation of chromosomes
		CHOPseq.MetaCorr(r[!is.na(r)], len[!is.na(r)], method='weighted.mean'); 		
	}
	#########################################

corr<-sapply(cov, function(cov, len) getCorr(cov, len), len=effective.len); 

individual<-data.frame(Score=scores, Count=score.count, Corr=corr); # return object
out<-list(effective.size=effective.len, individual=individual);

if (verbose) print('Finished getting individual mapq scores');

############################ run stepwise process
if (stepwise) {
if (verbose) print('Running stepwise process...');

if (toupper(stepwise.order)=='BEST') i<-3
else if (toupper(stepwise.order)=='LARGEST') i<-2
else i<-1;
s<-rownames(individual[order(individual[, i], decreasing=T), ]); 

# initiate stepwise process
sw.cov<-cov[[s[1]]]; # list of 2 coverage objects for both strands
sw.score<-s[1];
sw.count<-individual[s[1], 2];
sw.corr<-individual[s[1], 3];

remain<-s[-1];

    ####################### pmax 2 coverage objects
    pmaxCov<-function(cov1, cov2) lapply(1:length(cov1), function(i, cov1, cov2) pmax(cov1[[i]], cov2[[i]]), cov1=cov1, cov2=cov2);
    ###############################################
    
###########################
while(length(remain)>0) { # while there are still unused scores
if (stepwise.next.best) { # pick the next score that improves the strand correlation the most, very time comsuming
r.add<-sapply(remain, function(score, cov, sw.cov, len) 
	getCorr(list(pmaxCov(sw.cov[[1]], cov[[score]][[1]]), pmaxCov(sw.cov[[2]], cov[[score]][[2]])), len), 
cov=cov[remain], sw.cov=sw.cov, len=effective.len); 
next.score<-remain[which(r.add==max(r.add))[1]];
print(paste('Next score added:', next.score)); 
}
else next.score<-remain[1];

sw.cov[[1]]<-pmaxCov(sw.cov[[1]], cov[[next.score]][[1]]); 
sw.cov[[2]]<-pmaxCov(sw.cov[[2]], cov[[next.score]][[2]]);
sw.score<-c(sw.score, next.score);
sw.count<-c(sw.count, sum(individual[sw.score,2]));
sw.corr<-c(sw.corr, getCorr(sw.cov, effective.len));
if (verbose) print(paste('Added score=', next.score, ', corr=', sw.corr[length(sw.corr)], sep='')); 

remain<-remain[remain!=next.score];
} 
########################### end of while

stepwise<-data.frame(Score_Ordered=as.numeric(sw.score), Count_Total=sw.count, Corr=sw.corr);
out<-append(out, list(stepwise=list(results=stepwise, order.by=stepwise.order, next.best=stepwise.next.best)));  
} ########################## end of run stepwise process

out;
}
