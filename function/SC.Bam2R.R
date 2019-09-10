# load aligned read in BAM file, create R objects and metadata files
SC.Bam2R<-function(from, to, mapped, unmapped=NA, meta=list(), chr=NA, fields=NA) {
wd<-getwd();

if (file.exists(to)) setwd(to)
else {
dir.create(to);
setwd(to);
}

mapped<-paste(from, mapped, sep='/');
sub('//', '/', mapped)->mapped;
if (file.exists(mapped)) cat('Reading from file', mapped, '\n')
else stop('Bam file not exists');

if (identical(chr, NA)) chr<-names(scanBamHeader(mapped)[[1]][[1]])
else chr<-intersect(chr, names(scanBamHeader(mapped)[[1]][[1]]));

if (identical(chr, NA)) fields<-SC.SamBasicFields()
else fields<-union(SC.SamBasicFields(), fields);

fields<-intersect(fields, scanBamWhat());

tags<-SC.BamScannerChr(mapped, chr=chr, fields=fields, by.chr=TRUE, verbose=TRUE);

if (!is.na(unmapped)) {
unmapped<-paste(from, unmapped, sep='/');
sub('//', '/', unmapped)->unmapped;
if (file.exists(unmapped)) {
qname<-scanBam(unmapped, param=ScanBamParam(what='qname'))[[1]][[1]];
save(qname, file='unmapped.rdata', compress=TRUE);
}
}

chr<-unlist(lapply(tags, function(x) as.vector(x[,'rname'])))
strand<-unlist(lapply(tags, function(x) as.vector(x[,'strand'])));
pos<-unlist(lapply(tags, function(x) x[['pos']]));
width<-unlist(lapply(tags, function(x) x[['qwidth']]));

gr<-GRanges(seqnames=chr, strand=strand, ranges=IRanges(start=pos, width=width));
seqlen<-scanBamHeader(mapped)[[1]][[1]][seqlevels(gr)];
len<-seqlen[as.vector(seqnames(gr))];
gr<-gr[len>=end(gr)];
seqlengths(gr)<-seqlen;
metadata(gr)<-meta;

extra<-setdiff(fields, SC.SamBasicFields());

for (i in 1:length(extra)) {
assign(extra[i], unlist(lapply(tags, function(x) as.vector(x[, extra[i]])), use.names=FALSE), );
save(list=extra[i], file=paste('all_', extra[i], '.rdata', sep=''), compress=TRUE);
}

save(gr, file='all_aligned.rdata', compress=TRUE);

setwd(wd);
gr;
}
