# Summary read counts around given regions
SC.ReadCount<-function(query, subject, ws=0, no.dup=FALSE, adjust.strand=TRUE) {
# query, subject	GRanges object, input to findOverlaps function
# ws				Extension on both sides of the query regions
# no.dup			If TRUE, remove duplicated reads
# adjust.strand		Whether to determine relative position based on strand of query regions

subject<-resize(subject, 1);
query<-resize(query, 1);
#loc0<-start(query); 

str0<-as.vector(strand(query));
str1<-as.vector(strand(subject));
strand(query)<-strand(subject)<-'*';

#if (ws[1]>0) {
#start(query)<-start(query)-ws[1];
#end(query)<-end(query)+ws[1];
#}

olap<-as.matrix(findOverlaps(query, subject, maxgap=ws));
#loc1<-start(subject)[olap[,2]];

adj<-rep(1, nrow(olap));
if (adjust.strand) {
adj[str1[olap[,2]]=='-']<--1;
}
dist<-adj*(start(query)[olap[,1]]-start(subject)[olap[,2]]);

s<-split(dist, olap[,1]);
if (no.dup) s<-lapply(s, unique);

# read count at each query region
sum.row<-rep(0, length(query));
sum.row[as.numeric(names(s))]<-sapply(s, length);

# read count at each position relative to query regions
sm<-table(unlist(s));
sum.col<-rep(0, 1+2*ws);
names(sum.col)<-c((-1*ws):ws);
sum.col[names(sm)]<-sm;

list(row=sum.row, col=sum.col, all=s);
}
