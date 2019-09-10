# Retrieve the base sequences from a genome
SC.RetrieveSeq<-function(gr, genome, before=0, after=before, fix='center') {
# gr            GRanges, the regions around which the sequences to be retrieved
# genome        BSgenome, the sequences in a genome
# before, after Integer, number of bases to expand from the fixed points
# fix           String, 'start', 'end', 'both', or 'center', define the regions to retrieve

if (!(fix %in% c('start', 'end', 'both'))) fix<-'center';
before<-max(0, before);
after<-max(0, after);

# memorize the order
if (is.null(names(gr))) names(gr)<-1:length(gr)
ind<-1:length(gr)
names(ind)<-names(gr)

getSeq<-function(gr, chr, before, after, fix) {
len<-length(chr)
if (fix=='both') {
loc1<-start(gr)-before;
loc2<-end(gr)+after;
seq<-sapply(1:length(loc1), function(i, loc1, loc2, chr) as.character(subseq(chr, loc1[i], loc2[i])), loc1=pmax(1, loc1), loc2=pmin(len, loc2), chr=chr);
}
else {
if (fix=='start') loc<-start(gr)
else if (fix=='end') loc<-end(gr)
else loc<-start(gr)/2+end(gr)/2;
seq<-sapply(1:length(loc), function(i, loc1, loc2, chr) as.character(subseq(chr, loc1[i], loc2[i])), loc1=pmax(1, loc-before), loc2=pmin(len, loc+after), chr=chr);
}
names(seq)<-names(gr);
seq
}

chr<-unique(as.character(seqnames(gr)))
seq<-lapply(chr, function(chr, gr, g, bf, af) getSeq(gr[as.character(seqnames(gr))==chr], g[[chr]], bf, af, fix), gr=gr, g=genome, bf=before, af=after)
seq<-unlist(seq)
seq<-seq[order(ind[names(seq)])]; # re-order to match the order of input

cat('Retrieved ', length(seq), ' sequences', '\n')
seq
}