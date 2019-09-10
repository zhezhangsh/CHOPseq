SAM.BASIC.FIELDS<-function() c('rname', 'strand', 'pos', 'qwidth');

# A wrapper of the Rsamtools:bamScanner function to read in the fields needed for GRanges object from a BAM file
# faster and take less disk space
SC.BamScannerSlim<-function(bamFile, verbose=FALSE) {

# bamFile	character, path and name of bam file

library(Rsamtools);
library(GenomicRanges);

fields<-SAM.BASIC.FIELDS();
len<-scanBamHeader(file)[[1]][[1]];
#chr<-names(len);

####################################################################
param<-ScanBamParam(what=fields)
taken<-system.time(tags<-scanBam(bamFile, param=param)[[1]])[3];
####################################################################
if (verbose) cat('Read in', length(tags[[1]]), 'reads in', taken, 'seconds\n');

gr<-GRanges(seqnames=as.vector(tags$rname), strand=as.vector(tags$strand), ranges=IRanges(start=tags$pos, width=tags$qwidth));
#seqlengths(gr)<-len;


gr;
}
