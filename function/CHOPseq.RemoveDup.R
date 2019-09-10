# remove duplicate tags: tags on the same chromosome and strand and have the same start position
CHOPseq.RemoveDup<-function(gr, same.str=T) {
# gr		GRanges, tag collection to be filtered, tags are sorted by chromosome, then by start position, then by strand
# same.str	boolean, if true, only tags on the same strand are considered duplicates 
library(chipseq);

chr<-as.character(seqnames(gr)); 
if (same.str) str<-as.character(strand(gr));
start<-start(gr);

t<-chr[-1]!=chr[-length(chr)]; # if t is TRUE, a tag is not a duplicates of its previous one
if (same.str) t[str[-1]!=str[-length(str)]]<-T;
t[start[-1]!=start[-length(start)]]<-T;
t<-c(TRUE, t); # the first tags will always be kept

gr[t];
}