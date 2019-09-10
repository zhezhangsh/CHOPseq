# Perform a fast t test on a batch of variables, only get t statistics and p values
SC.FastT<-function(d, ind1, ind2, equal.var=FALSE) {

if (class(d)!='matrix') d<-as.matrix(d);

ind1<-intersect(1:ncol(d), ind1);
ind2<-intersect(1:ncol(d), ind2);

n1<-length(ind1);
n2<-length(ind2);

if (n1<=1|n2<=1) stop('Not enough observations');

d1<-d[, ind1];
d2<-d[, ind2];

m1<-rowMeans(d1);
m2<-rowMeans(d2);
diff<-m2-m1;

if (equal.var) {
ss1<-rowSums((d1-m1)^2)/(n1-1);
ss2<-rowSums((d2-m2)^2)/(n2-1);

t<-diff/(sqrt((rowSums((d1-m1)^2)+rowSums((d2-m2)^2))/(n1+n2-2))*sqrt(1/n1+1/n2));
df<-n1+n2-2;
}
else {
ss1<-rowSums((d1-m1)^2)/(n1-1);
ss2<-rowSums((d2-m2)^2)/(n2-1);
s<-sqrt(ss1/n1+ss2/n2);
t<-diff/s;
df<-(ss1/n1+ss2/n2)^2/((ss1/n1)^2/(n1-1)+(ss2/n2)^2/(n2-1));
}

p<-rep(0.5, length(t));
p[t>0]<-1-pt(t[t>0], df=df[t>0]);
p[t<0]<-pt(t[t<0], df=df[t<0]);
p<-p*2;

list(diff=diff, t=t, p=p);
}
