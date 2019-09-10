# Compare 2 sets of counts to generate a set of statistics
# For ChIPseq data, a common analysis is to compare the counts of reads in certain regions (such as promoters) of a pair of samples
SC.ComparePairOfCounts<-function(C0, C1, anno, names=c('C0', 'C1'), N0=0, N1=0, cutoff.p=1, cutoff.l2r=0, cutoff.c=0, num.top=nrow(anno), order.by=c('l2r', 'c', 'p'), p.method=c('proportion', 'poisson')) {
# C0, C1        Vectors of counts with the same length
# anno          data.frame, annotation of C0 and C1
# N0, N1        Total counts in both samples (number of all aligned reads for ChIPseq), used to normalize C0, C1; do normalization only if both values are greater than 0
# cutoff        Cutoff values to select entries of interest, respectively p value, log2-ratio and counts
# p.method      Statistical method used to generate p values
# num.top       Maximum number of top entries to select
# order.by        If num.top is greater than 0, how to select the top entries, by p value, count, or log2-ratio
# DAVID         Whether to run a DAVID analysis with the top entries

if (num.top>0) top<-order.by[1] else top<-'';
if (!(top %in% order.by)) top<-'';

id<-intersect(names(C0), names(C1));
id<-intersect(id, rownames(anno));
C0<-C0[id];
C1<-C1[id];
anno<-anno[id,];

if (length(id)<=1) {
    cat('No common IDs, return NA\n');
    NA;
} else {
results<-list();

l2r<-log2(pmax(C1, 1)/pmax(C0, 1));

# get p value using a given statistical test
if (N0>0 & N1>0 & p.method[1]=='proportion') p.method=='proportion' else p.method='poisson';
if (p.method[1]=='proportion') {
p<-apply(cbind(C0, C1), 1, function(c) prop.test(c, c(N0, N1))$p.value[[1]]);
}
else {
if (N0>0 & N1>0) r<-N1/N0 else r<-1;
p<-apply(cbind(C0, C1), 1, function(c) poisson.test(c, r=r)$p.value[[1]]);
}
p[C0==0&C1==0]<-1;
q<-p.adjust(p, method='fdr');

d<-data.frame(C0, C1, l2r, p, q);
names(d)[1:2]<-paste('Count_', names, sep='');
names(d)[3]<-'Log2-ratio';
names(d)[4]<-paste('P-', p.method, sep='');
names(d)[5]<-'FDR';

results[['Complete List']]<-cbind(anno, d);

up<-d[p<=cutoff.p&l2r>=abs(cutoff.l2r)&C1>=abs(cutoff.c), ];
dn<-d[p<=cutoff.p&l2r<=-1*abs(cutoff.l2r)&C0>=abs(cutoff.c), ];

if (order.by[1]=='p') {
up<-up[order(up[, 4]), ][1:min(nrow(up), num.top),];
dn<-dn[order(dn[, 4]), ][1:min(nrow(dn), num.top),];
} else if (order.by[1]=='c') {
up<-up[order(up[, 2]), ][nrow(up):(nrow(up)-min(nrow(up), num.top)+1),];
dn<-dn[order(dn[, 1]), ][nrow(dn):(nrow(dn)-min(nrow(dn), num.top)+1),];
} else {
up<-up[order(up[, 3]), ][nrow(up):(nrow(up)-min(nrow(up), num.top)+1),];
dn<-dn[order(dn[, 3]), ][1:min(nrow(dn), num.top),];
}

results[[paste('Higher in', names[1])]]<-cbind(anno[rownames(dn), ], dn);
results[[paste('Higher in', names[2])]]<-cbind(anno[rownames(up), ], up);

results;
} # end of else
}