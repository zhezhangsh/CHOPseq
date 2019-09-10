# Parse an entry of RefSeq gene downloaded from UCSC
CHOPseq.ParseRefGene<-function(gene) {
gene<-as.list(gene); 
id<-as.vector(gene[['ID']]);
symbol<-as.vector(gene[['Symbol']]);
chr<-as.vector(gene[['Chromosome']]);
strand<-as.vector(gene[['Strand']]);
if (strand==1) strand='+' else if (strand==-1) strand='-';

if (strand=='+'|strand=='-') { # strand must be '+' or '-'
# transcript start/end
start0<-1+as.numeric(gene[['TranscriptStart']]);
end0<-as.numeric(gene[['TranscriptEnd']]);
# coding sequence start/end
start1<-1+as.numeric(gene[['CdsStart']]);
end1<-as.numeric(gene[['CdsEnd']]);

# promoter and UTRs
if (strand=='+') 
regions<-list(c(start0-10000, start0-1), c(start0-1000, start0-1), c(start0, start1-1), c(end1+1, end0))
else 
regions<-list(c(end0+1, end0+10000), c(end0+1, end0+1000), c(end1+1, end0), c(start0, start1-1));
x<-do.call('rbind', regions);
x<-data.frame(Region=c('10K-Promoter', '1K-Promoter', '5-UTR', '3-UTR'), Start=x[,1], End=x[,2]);
# non-coding genes have no UTRs
if (start1==1+end1) x<-x[-(3:4),]

# exons
starts<-as.vector(gene[['ExonStarts']]);
ends<-as.vector(gene[['ExonEnds']]);
starts<-as.numeric(strsplit(starts, ',')[[1]])+1;
ends<-as.numeric(strsplit(ends, ',')[[1]]);
exon<-data.frame(Region=rep('Exon', length(ends)), Start=starts, End=ends);
x<-rbind(x, exon);

# if there is any intron
if (nrow(exon)>1) {
intron<-data.frame(Region=rep('Intron', length(ends)-1), Start=ends[-nrow(exon)]+1, End=starts[-1]-1);
x<-rbind(x, intron);
}

# coding regions
coding<-x[x[, 'Region']=='Exon',];
coding[coding[, 'Start']<start1, 'Start']<-start1;
coding[coding[, 'End']>end1, 'End']<-end1;
coding[, 'Region']<-'Coding';
x<-rbind(x, coding);

x<-x[x[, 'Start']<=x[, 'End'], ];
n<-nrow(x);
x<-data.frame(ID=rep(id, n), Symbol=rep(symbol, n), Chromosome=rep(chr, n), Strand=rep(strand,n), x); 

}
else x<-NA;

x;
}
