# Plot statistics of exons in a transcript
SC.PlotExons<-function(height, start, end, col=rep('grey', length(height)), main='', ylab='', xlab='Block size, Log10(bases)', legend=FALSE) {

# exon width
width<-abs(end-start)+1;
width[width<1]<-1;
width<-log10(width);

# distance betwee exons
gap<-c(1, abs(start[-1]-end[-length(end)]));
gap[gap<1]<-1;
gap<-log10(gap);

x1<-sapply(1:length(gap), function(i, gap, width) sum(gap[1:i])+sum(width[1:i]), gap=gap, width=width);
x0<-x1-width;

xlim<-c(0, sum(gap+width));
min<-min(height);
max<-max(height);
range<-max-min;
if (legend) ylim<-c(max-1.05*range, min+1.2*range) else ylim<-c(max-1.05*range, min+1.1*range)

if (ylim[1]>0) ylim[1]<-0;
if (ylim[2]<0) ylim[2]<-0;

plot(0, type='n', xlim=xlim, ylim=ylim, main=main, ylab=ylab, col.axis=4, xaxt='n', xaxs='i'); 
box();
abline(h=0, col=8);
rect(x0, 0, x1, height, col=col, border=4);
title(xlab=xlab, line=1);
#axis(2, col=4, col.ticks=4);


#plot legend
if (legend) {
yaxp<-par()$yaxp;
yloc<-min+range*c(1.1, 1.15, 1.2);
xloc<-2:4;
rect(xloc, yloc, 6, yloc-.04*range, col=col);
text(6, yloc-0.02*range, pos=4, labels=c('10,000 bases', '1,000 bases', '100 bases'));
}

cbind(x0, x1);
}
