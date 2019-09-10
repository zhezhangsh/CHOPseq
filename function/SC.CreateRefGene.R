# Create an object that includes different types of features of RefSeq genes
# The RefSeq gene information will be downloaded from UCSC
SC.CreateRefGene<-function(genome.name) {
library(rtracklayer);

session<-browserSession();
genome(session)<-genome.name;
query<-ucscTableQuery(session, 'refGene');
ref<-getTable(query);

if (toupper(genome.name)=='HG19') {
library(BSgenome.Hsapiens.UCSC.hg19);
gn<-Hsapiens;
CHR<-names(gn)[1:25];
}
else if (toupper(genome.name)=='HG18') {
library(BSgenome.Hsapiens.UCSC.hg18);
gn<-Hsapiens;
CHR<-names(gn)[1:25];
}
else if (toupper(genome.name)=='MM9') {
library(BSgenome.Mmusculus.UCSC.mm9);
gn<-Mmusculus;
CHR<-names(gn)[1:22];
}
else if (toupper(genome.name)=='MM8') {
library(BSgenome.Mmusculus.UCSC.mm8);
gn<-Mmusculus;
CHR<-names(gn)[1:22];
}
else stop('Error: Unknown genome name', genome.name);

chr.len<-seqlengths(gn)[CHR];

ref<-ref[as.vector(ref[, 'chrom']) %in% CHR, ];
tx<-GRanges(seqnames=as.vector(ref[, 'chrom']), ranges=IRanges(ref[, 'txStart'], ref[, 'txEnd']), strand=as.vector(ref[, 'strand']));
names(tx)<-1:length(tx);
seqlengths(tx)<-seqlengths(gn)[seqlevels(tx)];
fields<-c('name', 'name2', 'score', 'exonCount', 'cdsStart', 'cdsEnd', 'cdsStartStat', 'cdsEndStat');
fields<-fields[fields %in% colnames(ref)];
meta<-ref[, fields];
names(meta)[1:2]<-c('ID', 'Symbol');
for (i in 1:ncol(meta)) if (class(meta[[i]])=='factor') meta[[i]]<-as.vector(meta[[i]]);
elementMetadata(tx)<-meta;

refGene<-list(Transcripts=tx, TSS=resize(tx, 1));

cat('Transcripts added\n');

#Unique TSS
x<-paste(as.vector(ref[, 'chrom']), ref[, 'txStart'], ref[, 'txEnd'], as.vector(ref[, 'strand']), sep='_');
id<-as.vector(meta[,1]);
syml<-as.vector(meta[,2]);
id<-sapply(split(as.vector(meta[[1]]), x), function(a) a[1]);
symb<-sapply(split(as.vector(meta[[2]]),x),function(a) a[1]); 
y<-do.call('rbind', strsplit(names(id), '_'));
tss<-GRanges(seqnames=as.vector(y[,1]),ranges=IRanges(as.numeric(as.vector(y[,2])),as.numeric(as.vector(y[,3]))),strand=as.vector(y[, 4]));
seqlengths(tss)<-chr.len[seqlevels(tss)];
elementMetadata(tss)$ID<-as.vector(id);
elementMetadata(tss)$Symbol<-as.vector(symb);
tss<-tss[order(as.vector(seqnames(tss)), start(tss))];
refGene$TSS.Unique<-tss;

cat('TSS added\n');

# Masks
masks<-lapply(CHR, function(c, gn) as.data.frame(masks(gn[[c]])) , gn=gn);
c<-rep(CHR, sapply(masks, nrow));
masks<-do.call('rbind', masks);
refGene$Masks<-GRanges(seqnames=c, ranges=IRanges(masks[,2], masks[,3]));
elementMetadata(refGene$Masks)$Type<-as.vector(masks[,1]);
seqlengths(refGene$Masks)<-seqlengths(gn)[seqlevels(refGene$Masks)];

cat('Masks added\n');

###########################################
# Sub-gene regions
#Exon-UTRs
stt<-as.vector(ref[, 'exonStarts']);
end<-as.vector(ref[, 'exonEnds']);
stt<-lapply(strsplit(stt, ','), as.numeric);
end<-lapply(strsplit(end, ','), as.numeric);
cds.stt<-ref[, 'cdsStart'];
cds.end<-ref[, 'cdsEnd'];
str<-as.vector(ref[, 'strand']);
regions<-lapply(1:length(str), function(i) SC.ParseSubRegions(stt[[i]], end[[i]], cds.stt[i], cds.end[i], str[i]));
n<-sapply(regions, nrow);
regions<-do.call('rbind', regions);
mt<-data.frame(Type=regions[,3], ID=rep(as.vector(ref[,'name']), n), Symbol=rep(as.vector(ref[,'name2']), n));
regions<-GRanges(seqnames=rep(as.vector(ref[, 'chrom']), n), ranges=IRanges(regions[,1], regions[,2]), strand=rep(as.vector(ref[,'strand']),n));

# promoters
prm<-resize(refGene$TSS, 1);
prm<-resize(prm, 2, fix='end');
prm<-resize(prm, 1);
prm1k<-as.data.frame(resize(prm, 1000, fix='end'));
prm10k<-as.data.frame(resize(prm, 10000, fix='end'));
prm1k$Type<-rep('1K-Promoter', nrow(prm1k));
prm10k$Type<-rep('10K-Promoter', nrow(prm10k));
prm<-rbind(prm1k, prm10k);
regions<-c(regions, GRanges(seqnames=as.vector(prm[,1]), ranges=IRanges(prm[,2], prm[,3]), strand=as.vector(prm[,5])));
mt<-rbind(mt, prm[, c('Type', 'ID', 'Symbol')]);

# 1k downstream
dn1k<-resize(refGene$TSS, 1000);
elementMetadata(dn1k)<-NULL;
regions<-c(regions, dn1k);
mt<-rbind(mt, data.frame(Type=rep('1K-TSS-Downstream',length(dn1k)),ID=elementMetadata(refGene$TSS)$ID, Symbol=elementMetadata(refGene$TSS)$Symbol));

#Assembly gaps
m<-refGene$Masks;
gaps<-as.data.frame(m[elementMetadata(m)[[1]]=='AGAPS']);
regions<-c(regions, GRanges(seqnames=as.vector(gaps[,1]), ranges=IRanges(gaps[,2], gaps[,3]), strand=as.vector(gaps[,5])));
mt<-rbind(mt, data.frame(Type=rep('Gaps', nrow(gaps)), ID=rep('', nrow(gaps)), Symbol=rep('', nrow(gaps))));

# Intergenic regions
c<-coverage(regions);
v<-cbind(unlist(lapply(c, runValue)), unlist(lapply(c, function(c) start(c))), unlist(lapply(c, function(c) end(c))));
chr<-rep(names(c), sapply(c, nrun))[v[,1]==0];
v<-v[v[,1]==0,];
regions<-c(regions, GRanges(seqnames=chr, ranges=IRanges(v[,2], v[,3])));
mt<-rbind(mt, data.frame(Type=rep('Intergenic', nrow(v)), ID=rep('', nrow(v)), Symbol=rep('', nrow(v))));

names(mt)[1]<-'Region';
elementMetadata(regions)<-mt;
refGene$Regions<-regions;
cat('Sub-gene regions added\n');
############################################################################################

refGene;
}


# Split a refSeq genes into sub-gene regions
SC.ParseSubRegions<-function(exon.stt, exon.end, cds.stt, cds.end, strand) {
# exon.stt, exon.end	Numeric vectors, Start and end positions of all exons, ordered from low to high
# cds.stt, cds.end	Numeric, Start and end position of coding region
# strand		Value of -1, 1, '+', or '-'

# the lowest position
tss<-exon.stt[1]-1;

exon.stt<-exon.stt-tss;
exon.end<-exon.end-tss;
cds.stt<-cds.stt-tss;
cds.end<-cds.end-tss;
 
c<-coverage(GRanges(seqnames='1', ranges=IRanges(exon.stt, exon.end)))[[1]];
c[cds.stt:cds.end]<--1*c[cds.stt:cds.end];

stt<-start(c)+tss;
end<-end(c)+tss;
v<-runValue(c)+2;

type<-c('Coding', 'Intron', 'UTR');
type<-type[v];

ind<-which(type=='Coding');
ind0<-which(type=='UTR');

# specify UTR type
if (strand==1 | strand=='+') {
type[ind0[ind0<ind[1]]]<-'5-UTR';
type[ind0[ind0>ind[length(ind)]]]<-'3-UTR';
}
else if (strand==-1 | strand=='-') {
type[ind0[ind0<ind[1]]]<-'3-UTR';
type[ind0[ind0>ind[length(ind)]]]<-'5-UTR';
}

regions<-data.frame(Start=stt, End=end, Type=type);
exon<-data.frame(Start=exon.stt+tss, End=exon.end+tss, Type=rep('Exon', length(exon.stt)));

regions<-rbind(regions, exon);

# no coding region, no UTRs
if ((cds.end-cds.stt)<=0) {
regions<-regions[regions[,3]!='Coding',];
regions<-regions[regions[,3]!='3-UTR',];
regions<-regions[regions[,3]!='5-UTR',];
}

regions<-regions[(regions[,2]-regions[,1])>0,]; # remove 0-base long regions

regions;
}
