# Match a PWM to a set of chromosomes and format the output as a GRanges
CHOPseq.MatchPWM<-function(pwm, genome, chr, min.score='', min.dist=0, strand=0) {
# pwm		Data matrix, a position weight matrix with 4 rows named 'A', 'C', 'G', 'T'
# genome	BSgenome, the whole sequence information of genome
# chr       String vector, the name of the chromosomes to search the PWM
# min.score	String, used by matchPWM(Biostrings) function
# min.dist	Non-negative Integer, the minimum distance allow between 2 hits 
# strand	Integer, if 0, both strand; if >0, forward strand only; if <0, reverse strand only
chr<-intersect(names(genome), chr); 

if (length(chr)>0) {
hits<-lapply(chr, function(chr, pwm, s, g) CHOPseq.MatchPWMChr(pwm, g[[chr]], min.score=s), pwm=pwm, s=min.score, g=genome);
chr<-rep(chr, sapply(hits, nrow));
hits<-do.call('rbind', hits);
s<-c('-', '+');
strand<-s[(hits[,3]+3)/2];
GRanges(seqnames=chr, ranges=IRanges(start=hits[,1], end=hits[,2]), strand=strand, seq=as.vector(hits[,4]));
}
else NA;
}


# Match a PWM to a chromosome and format the output as a sorted data frame
CHOPseq.MatchPWMChr<-function(pwm, seq, min.score='', min.dist=0, strand=0) {
# pwm		Data matrix, a position weight matrix with 4 rows named 'A', 'C', 'G', 'T'
# seq		Biostrings, the whole sequence information of a chromosome, available from a BSgenome
# min.score	String, used by matchPWM(Biostrings) function
# min.dist	Non-negative Integer, the minimum distance allow between 2 hits 
# strand	Integer, if 0, both strand; if >0, forward strand only; if <0, reverse strand only
library(Biostrings);

	# the working horse, match the PWM to chromosome and format the outputs
	MatchAndFormat<-function(pwm, seq, min.score, str) {
	matched<-matchPWM(pwm, seq, min.score=min.score);
	data.frame(start=start(matched), end=end(matched), strand=rep(str, length(matched)), seq=as.character(matched));
	} 

if (strand>=0) forward<-MatchAndFormat(pwm, seq, min.score, 1);
if (strand<=0) reverse<-MatchAndFormat(reverseComplement(pwm), seq, min.score, -1); 

if (strand>0) out<-forward # forward strand only
else if (strand<0) out<-reverse # reverse strand only
else { # if both strand
	out<-rbind(forward, reverse); 
	out<-out[order(out[,1]), ];
};

# Perform this step to avoid having too many matched locations in a small region
if (min.dist>0) { 
dist<-c(min.dist, out[-1, 1]-out[-nrow(out), 1]); 
out<-out[dist>=min.dist, ]; # if the distance between 2 locations is smaller than min.dist, only keep the first one
# todo: by strand?
}

out; 
} # end of function
# Return, a data frame with 4 columns: start, end, strand (-1/1), and DNA sequence
