# Count number of aligned locations of reads with different length
# use snow package to distribute tasks
SC.ClusterAlignCount<-function(seq, reference, NC=8, mismatch=0, verbose=TRUE) {
#seq		String vector, the read sequence to be aligned
#reference	Name of the reference genome
#NC		Number of Snow clusters
#mismatch	The number of mismatch bases
if (length(seq)==0) {
cat('There are no reads to be aligned\n');
NA;
}
else {
library(BSgenome);

if (toupper(reference)=='HG19') {
library(BSgenome.Hsapiens.UCSC.hg19);
n.chr<-25;
CHR<-names(Hsapiens)[1:n.chr];
genome<-Hsapiens;
}
else if (toupper(reference)=='MM9') {
library(BSgenome.Mmusculus.UCSC.mm9);
n.chr<-22;
CHR<-names(Mmusculus)[1:n.chr];
genome<-Mmusculus;
}
else if (toupper(reference)=='TEST') {
library(BSgenome.Hsapiens.UCSC.hg19);
n.chr<-2
CHR<-names(Hsapiens)[c(9, 21)];
genome<-Hsapiens;
}
else {
stop(paste('Unknown genome ', genome.name, sep=''));
}

if (class(seq)!='DNAStringSet') seq<-DNAStringSet(seq);

if (verbose) cat('Splitting reads based their width\n');
width<-width(seq);
#ind<-split(1:length(seq), width(seq));
#ind<-ind[length(ind):1];
seq<-split(seq, width); 

# create Clusters
library(snow);
cl<-makeCluster(NC, type='SOCK');
if (verbose) cat('Aligning through threading\n');
aligned<-clusterApplyLB(cl, seq, vcountPDict, genome, exclude=setdiff(names(genome), CHR), 
	max.mismatch=mismatch, min.mismatch=mismatch);
stopCluster(cl);

if (verbose) cat('Alignment done\n');
aligned<-lapply(aligned, function(x) x[[4]]);
count<-lapply(aligned, function(x) rowSums(matrix(as.vector(x), nr=length(x)/2/length(CHR))));
count<-unsplit(count, width);
#count<-unlist(count, use.names=FALSE);
#ind<-unlist(ind, use.names=FALSE);
#count[ind]<-count;
count;
} # end of else
}

