# From coverage data, set depth of islands  whose max depth lower than a given cutoff to zero
SC.RleRemoveBG.R<-function(rle, min.height=5, lower=1) {
# rle	RleList, depth info.
# min.height	Minimum peak height to keep
# lower	Lower cutoff for slicing islands

lapply(rle, function(c) {
slice<-slice(c, lower=lower);
if (length(slice)==0) {
c;
}
else {
height<-viewMaxs(slice);
seqselect(c, start(slice)[height<min.height], end(slice)[height<min.height])<-0;
c;
}
})
}
