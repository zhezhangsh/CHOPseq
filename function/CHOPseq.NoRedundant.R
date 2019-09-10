# Remove rows of a data frame or matrix have the same values in all given columns
CHOPseq.NoRedundant<-function(x, columns=1:ncol(x)) {
c<-intersect(1:ncol(x), columns);
y<-as.data.frame(x[, c]);
z<-lapply(y, as.vector);
a<-do.call('paste', z);
b<-split(1:length(a), a);
d<-sapply(b, function(b) b[1]);
x[d,];
}
