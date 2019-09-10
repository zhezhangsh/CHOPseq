# Perform a 2 group comparison with SAM method
CHOPbic.SAM.2Groups<-function(d, ind.g1=1:round(ncol(d)/2), ind.g2=(max(ind.g1)+1):ncol(d), groupnames=c('A', 'B'), paired=FALSE, nperm=100) {
# d		Data.frame or matrix
# ind.g1, ind.g2	Column indexes of the 2 groups
# groupnames		Names of the groups
# paired			Paired test?
# nperm				SAM parameter

library(samr);

if (paired) {
s<-c(-1*(1:length(ind.g1)), 1:length(ind.g2)); 
type<-"Two class paired"
}
else {
s<-rep(1:2, c(length(ind.g1), length(ind.g2)));
type<-"Two class unpaired";
}

# run SAM and format outputs
data<-list(x=cbind(d[, ind.g1], d[, ind.g2]), y=s, geneid=rownames(d), genenames=rownames(d), logged2=TRUE); 

capture.output(sam<-samr(data, resp.type=type, nperms=nperm));
p<-samr.pvalues.from.perms(sam$tt, sam$ttstar);
capture.output(delta<-samr.compute.delta.table(sam));
x<-samr.compute.siggenes.table(sam, 0, data, delta, all.genes=T); 
x<-rbind(x[[1]], x[[2]]);  
rownames(x)<-as.vector(x[,2]);  
x<-x[rownames(d), ]; 
fc<-as.numeric(x[,7]); 
q<-as.numeric(x[,8])/100; 
m1<-rowMeans(d[, ind.g1]); 
m2<-rowMeans(d[, ind.g2]);
nm<-paste(groupnames[1], groupnames[2], sep='-');
out<-cbind(m1, m2, m2-m1, fc, p, q);
colnames(out)<-c(paste('Mean_', groupnames[1], sep=''), paste('Mean_', groupnames[2], sep=''), nm, 'Fold_change', 'p_SAM', 'FDR_SAM');

out;
}