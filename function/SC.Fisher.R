# Run a fisher test
SC.Fisher<-function(s1, s2, u) {
# s1, s2    2 sets of values
# u         the universe


if (length(u)>=2) {
s1<-intersect(s1, u);
s2<-intersect(s2, u);
}

if (length(s1)==0&length(s2)==0) {
cat('Both sets are empty.\n');
NA;
} else {
s0<-intersect(s1, s2);
l0<-length(s0);
l1<-length(s1)-l0;
l2<-length(s2)-l0;
l<-length(u)-l0-l1-l2;
fisher.test(matrix(c(l0, l1, l2, l), nr=2));
} # end of else

}
