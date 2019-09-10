# Remove 'N' from sequence reads
SC.RemoveN<-function(seq) {
n<-grep('N', seq);
has.n<-length(n);

# remove end N
seq0<-seq[n];
while(has.n>0) {
cat(has.n, '\n');
first<-grep('^N', seq0);
last<-grep('N$', seq0);
seq0[first]<-substr(seq0[first], 2, nchar(seq0[first]));
seq0[last]<-substr(seq0[last], 1, nchar(seq0[last])-1);
if (has.n==length(grep('N', seq0))) break else has.n<-length(grep('N', seq0));
}

# if there is N in the middle of reads, use the longest segment without N
middle<-seq0[grep('N', seq0)];
if (length(middle)>0) {
middle<-strsplit(middle, 'N');
middle<-sapply(middle, function(x) {y<-nchar(x); x[y==max(y)][1];})
seq0[grep('N', seq0)]<-middle;
}
seq[n]<-seq0; # all Ns removed
seq;
}
