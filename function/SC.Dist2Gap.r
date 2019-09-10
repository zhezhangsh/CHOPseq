# get the nearest distance of reads to any assembly gaps
SC.DistToGap<-function(file, genome.name) {
gr<-eval(parse(text=load(file))); 

if (toupper(genome.name)=='HG19') refGene<-eval(parse(text=load('~/R/data/hg19/RefGene.RData')))
else if (toupper(genome.name)=='HG18') refGene<-eval(parse(text=load('~/R/data/hg18/RefGene.RData')))
else if (toupper(genome.name)=='MM9') refGene<-eval(parse(text=load('~/R/data/mm9/RefGene.RData')))
else if (toupper(genome.name)=='MM8') refGene<-eval(parse(text=load('~/R/data/mm8/RefGene.RData')))
else {}

masks<-refGene$Masks;
gaps<-masks[elementMetadata(masks)[[1]]=='AGAPS'];

dist<-rep(-1, length(gr));
c<-as.vector(seqnames(gr));
cg<-as.vector(seqnames(gaps));
sl<-seqlevels(gr);

for (j in 1:length(sl)) {
dist[c==sl[j]]<-SC.Dist2Nearest(gr[c==sl[j]], gaps[cg==sl[j]]);
}

elementMetadata(gr)$dist2gap<-dist+1;

gr;
}
