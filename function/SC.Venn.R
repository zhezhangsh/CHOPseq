# Draw a 2-set Venn diagram
SC.Venn<-function(s1, s2, names=c('Set1', 'Set2'),  universe=c(), fisher=TRUE) {

if (length(universe)>=2) {
s1<-intersect(s1, universe);
s2<-intersect(s2, universe);
}

quartz(w=8, h=6);

if (length(s1)==0&length(s2)==0) {
cat('Both sets are empty.\n');
NA;
} else {
out<-list();
par(mar=rep(0, 4), omi=rep(0, 4)); # no margins
plot(0, type='n', xlim=c(0, 8), ylim=c(0, 6), axes=F, bty='n', xaxs='i', yaxs='i', xlab='', ylab=''); # plot an empty space
lwd<-12;
bg<-'yellow1';
fg<-'seagreen3';
symbols(3, 3, circles=2, bg=bg, fg=fg, lwd=lwd, add=TRUE, inches=FALSE);
symbols(5, 3, circles=2, bg=bg, fg=fg, lwd=lwd, add=TRUE, inches=FALSE);
symbols(3, 3, circles=2, fg=fg, lwd=lwd, add=TRUE, inches=FALSE);

rect(0, 0, 8, 6, lwd=lwd, border='deepskyblue2');

s0<-intersect(s1, s2);
l0<-length(s0);
l1<-length(s1)-l0;
l2<-length(s2)-l0;
l<-length(universe)-l0-l1-l2;
out$counts<-c(l, l1, l2, l0);
cex=3;
col='darkred';
text(2, 3, labels=l1, col=col, cex=cex);
text(6, 3, labels=l2, col=col, cex=cex);
text(4, 3, labels=l0, col=col, cex=cex);

text(0.5, 3, labels=names[1], cex=cex-1, col='darkblue', srt=90);
text(7.5, 3, labels=names[2], cex=cex-1, col='darkblue', srt=270);

if (length(universe)>=2) {
text(4, .5, labels=l, cex=cex, col=col);
if (fisher) {
out$fisher<-fisher.test(matrix(out$counts, nr=2));
p<-out$fisher[[1]];
ci<-out$fisher[[2]];
or<-out$fisher[[3]];

if (p<0.0001) p<-format(p, scientific=T, digit=2) else p<-round(p, 4);
p<-gsub(' ', '', p);
line<-paste(
paste('Odds ratio=', round(or, 2), sep=''), '; ',
paste('95% C.I.=(', round(ci[1], 2), ', ', round(ci[2], 2), ')', sep=''), '; ',
paste('p=', p), sep='');
text(4, 5.5, labels=line, cex=1.75);

} # end of fisher test
}
} # end of else

}
