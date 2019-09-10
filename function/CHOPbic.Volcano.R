# Draw a volvano plot with given group difference and p values
CHOPbic.Volcano<-function(d, p, names=c('A', 'B'), screen=TRUE, pdf=FALSE, tiff=FALSE, vlines=0, hlines=1) {
# d         values of group difference on x-axis
# p         p values, the same length as d
# screen    plot to the screen?
# pdf       write pdf file?
# tiff      write tiff file?
# vline     position of vertical line(s) to draw
# hline     position of horizontal line(s) to draw

p[p>1]<-1;
p[p<1e-300]<-1e-300;
y<--1*log10(p);
y.lim=max(y)+0.1;
x.lim=max(abs(d));

size<-sqrt(abs(d)*y);
size<-2*size/max(size);

y.axis<-exp((0:round(y.lim))*log(.1));
y.axis<-sapply(1:length(y.axis), function(i, y) round(y[i], i-1), y=y.axis);


# write to pdf file
if (pdf) {
filename<-paste('Volcano, ', names[1], ' vs. ', names[2], '.pdf', sep='');
pdf(filename, w=10, h=7.5);

#########plotting
par(mar=c(5, 5, 2, 2));
plot(d, y, xlim=c(-1*x.lim, x.lim), ylim=c(0, y.lim), yaxs='i', yaxt='n', col='#333333', pch=19, cex=size, xlab='', ylab='');
if (!is.na(vlines)) abline(v=vlines, lty=3);
if (!is.na(hlines)&hlines<=1&hlines>=0) abline(h=-1*log10(hlines), lty=3);
axis(2, at=0:round(y.lim), labels=y.axis, las=2);
title(ylab='p value', line=3.2);
title(xlab=paste('Log2(', names[1], '/', names[2], ')', sep=''));
# end of plotting

capture.output(x<-lapply(dev.list(), function(x) dev.off(which=x)));
}

# write to tiff file
if (tiff) {
filename<-paste('Volcano, ', names[1], ' vs. ', names[2], '.tiff', sep='');
tiff(filename, w=1000, h=750);

#########plotting
par(mar=c(5, 5, 2, 2));
plot(d, y, xlim=c(-1*x.lim, x.lim), ylim=c(0, y.lim), yaxs='i', yaxt='n', col='#333333', pch=19, cex=size*1.25, xlab='', ylab='');
if (!is.na(vlines)) abline(v=vlines, lty=3);
if (!is.na(hlines)&hlines<=1&hlines>=0) abline(h=-1*log10(hlines), lty=3);
axis(2, at=0:round(y.lim), labels=y.axis, las=2);
title(ylab='p value', line=3.2, cex.lab=1.5);
title(xlab=paste('Log2(', names[1], '/', names[2], ')', sep=''), cex.lab=1.5);
# end of plotting

capture.output(x<-lapply(dev.list(), function(x) dev.off(which=x)));
}

# plot to screen
if (screen) {
#########plotting
par(mar=c(5, 5, 2, 2));
plot(d, y, xlim=c(-1*x.lim, x.lim), ylim=c(0, y.lim), yaxs='i', yaxt='n', col='#333333', pch=19, cex=size, xlab='', ylab='');
if (!is.na(vlines)) abline(v=vlines, lty=3);
if (!is.na(hlines)&hlines<=1&hlines>=0) abline(h=-1*log10(hlines), lty=3);
axis(2, at=0:round(y.lim), labels=y.axis, las=2);
title(ylab='p value', line=3.2);
title(xlab=paste('Log2(', names[1], '/', names[2], ')', sep=''));
# end of plotting
}

dev.cur();
}