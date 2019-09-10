# Perform a meta-analysis of correlation coefficients
SC.MetaCorr<-function(corr, n=rep(4, length(corr)), method='weighted.mean') {
# corr		Numeric vector, a set of correlation coefficients obtained from different tests
# n			Non-negative Integer vector, the sample size of different tests, same length as corr
# method	Character, the method used to summarized corr
#					"weighted.mean", default, mean of chromosomes weighted by effective chromosome sizes
#					"mean", simple mean, not recommended
#					"z", Fisher's r to z conversion. See www.statisticshell.com/meta.pdf or Field, A.P. (2001) Psychological Methods, 6(2), 161-180 for details

if (length(corr)!=length(n)) n<-rep(4, length(corr)); 

corr[abs(corr)>1]<-0; # coefficients range [-1, 1];
n[n<0]<- 0; # sample size >=0
n<-round(n); 

if (toupper(method)=='MEAN') R<-mean(corr)
else if (toupper(method)=='Z') {
Z<-weighted.mean(log((1+corr)/(1-corr))/2, w=n-3); 
R<-(exp(2*Z)-1)/(exp(2*Z)+1);
} 
else R<-weighted.mean(corr, w=n); 

R;
}


