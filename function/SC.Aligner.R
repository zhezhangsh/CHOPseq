SC.Aligner<-function(seq, reference, just.count=TRUE, min.mismatch=0, max.mismatch=0, max.reads=10, verbose=FALSE) {
#seq		String vector, the read sequence to be aligned
#reference	Name of the reference genome
#just.count	If TRUE, just count the number of best hits
#mismatch	The maximum number of mismatch base
library(BSgenome);

#if (toupper(reference)=='HG19') {
#library(BSgenome.Hsapiens.UCSC.hg19);
#n.chr<-25
#CHR<-names(Hsapiens)[1:n.chr];
#genome<-Hsapiens;
#}
#else {
#stop(paste('Unknown genome ', genome.name, sep=''));
#}

if (class(seq)!='DNAStringSet') seq<-DNAStringSet(seq);

min.mismatch<-max(0, min.mismatch);
max.mismatch<-max(min.mismatch, max.mismatch);

width<-width(seq);
min<-min(width);
max<-max(width);

# indexes of reads with given length
ind<-lapply(min:max, function(l, width) which(width==l), width=width);
names(ind)<-min:max;

# only support counting hits currently
if (!just.count) stop('just.count=FALSE is not supported yet')
else { # just count the number of hits
count<-data.frame(N.Mismatch=rep(0, length(seq)), N.Hits=rep(0, length(seq)));

for (i in 1:length(ind)) {
cat('Aligning', length(ind[[i]]), 'reads, length =', names(ind)[i], '\n');
if (length(ind[[i]])>0) count[ind[[i]],]<-SC.CountAlignFixed(seq[ind[[i]]], reference, min.mismatch, 
max.mismatch, max.reads, verbose);
#cat('Aligned reads with length =', names(ind)[i], '\n');
}

out<-count;
}
out;
}

SC.CountAlignFixed<-function(seq, reference, min.mismatch=0, max.mismatch=0, max.reads=10, verbose=FALSE) {
# seq	DNAStringSet, sequences to be aligned, all sequence should have the same length
# reference	Name of the reference genome
# mismatch	Maximum number of mismatches allowed
# max.reads	Maximum number of reads to be aligned altogether (in millions)
# return	The counts of hits with best match

if (toupper(reference)=='HG19') {
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
n.chr<-25
CHR<-names(Hsapiens)[1:n.chr];
genome<-Hsapiens;
}
else if (toupper(reference)=='TEST') {
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
n.chr<-2
CHR<-names(Hsapiens)[c(9, 21)];
genome<-Hsapiens;
}
else {
stop(paste('Unknown genome ', genome.name, sep=''));
}

if (class(seq)!='DNAStringSet') seq<-DNAStringSet(seq);
max.reads<-max.reads*1000000;

width<-width(seq);
if (max(width)!=min(width)) stop('Sequences have different length');

count<-data.frame(N.Mismatch=rep(0, length(seq)), N.Hits=rep(0, length(seq)));

# reads with ambiguous bases will be excluded
x<-alphabetFrequency(seq)[,'N'];
count[x>0, 1]<--1

#####################################################
for (i in min.mismatch:max.mismatch) {
ind.all<-which(count[,1]>=0&count[,2]<=0);
if (length(ind.all)==0) break else {
batches<-ceiling(length(ind.all)/max.reads);

for (b in 1:batches) {
ind<-ind.all[(b*max.reads-max.reads+1):min(length(ind.all), b*max.reads)];

cat('Aligning', length(ind), 'reads, mismatch =', i, '\n');
t<-system.time(x<-vcountPDict(seq[ind], genome, exclude=setdiff(names(genome), CHR), max.mismatch=i, min.mismatch=i)[,4]);
if (verbose) print(t);
c<-rowSums(matrix(as.vector(x), nr=length(ind)));
count[ind,1]<-i;
count[ind,2]<-c;
#cat('Aligned reads, mismatch =', i, '\n');
}}}

count;
}

