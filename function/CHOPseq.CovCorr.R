# calculated the pairwise correlation of coverage of a group of samples
CHOPseq.CovCorr<-function(cov, chr=NA) {
if (identical(NA, chr)) chr<-names(cov[[1]]);
cov<-lapply(cov, function(cov, chr) cov[chr], chr=chr);

es<-rowMeans(sapply(cov, CHOPseq.EffectiveSize));

r<-matrix(nr=length(cov), nc=length(cov));
rownames(r)<-colnames(r)<-names(cov);

for (i in 1:(length(cov)-1)) {
r[i,i]<-1;
for (j in (i+1):length(cov)) r[i, j]<-r[j, i]<-weighted.mean(cor(cov[[i]], cov[[j]]), w=es);
}
r[length(cov), length(cov)]<-1;
r;
}


# calculated the correlation of coverage of a pair of samples
CHOPseq.CovCorrPair<-function(cov1, cov2, chr=NA, es=NA) {
if (identical(NA, chr)) chr<-names(cov1);
if (identical(NA, es)) es<-CHOPseq.EffectiveSize(cov1+cov2)
else es<-es[chr];

cov1<-cov1[chr];
cov2<-cov2[chr];

r<-cor(cov1, cov2);

weighted.mean(r, w=es);
}
