# Evaluate uniqueness of sequences aligned by reads
SC.AlignmentRepetitive<-function(path, genome) {
#path   Full file names of a GRanges object saving the aligned reads
#genome Genome name, such as 'hg19'
  
gr<-eval(parse(text=load(path)));
seq<-getSeq(genome, gr);
seq<-SC.RemoveN(seq);
count<-SC.ClusterAlignCount(seq, genome@provider_version, mismatch=0);
elementMetadata(gr)$repetitive<-count;
#save(gr, file=path);

gr;
}
