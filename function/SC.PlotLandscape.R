# make a landscape plot of coverage depth on one or more chromosomes
SC.PlotLandscape<-function(depth, step=10000, width=step, plot=TRUE, pdf=FALSE, ylim='3SD', same.ylim=TRUE) {
# depth	a list of Rle objects, each represents a chromosome
# step	step size of scanning
# width	windows size to take the average depth
# plot plot landscape if TRUE
# pdf	write to pdf if TRUE
# ylim	up limit of y-axis; if '3SD', 3 SDs more than mean; if 'MAX', max value; 
# same.ylim	whether to use the same y-axis limit for each chromosome

	# calculate XY coordinates, X is the chromosome location and Y is mean depth on a chromosome
	getXY<-function(dep, step, width) {
		start<-seq(1, length(dep), step);
		end<-pmin(length(dep), start+width-1);
		list(X=start/2000000+end/2000000, Y=aggregate(dep, start=start, end=end, FUN=mean));
	}
		
XY<-lapply(depth, function(depth) getXY(depth, step, width));

if (plot) {
if (identical(toupper(ylim), '3SD')) ylim<-sapply(XY, function(x) mean(x[[2]])+3*sd(x[[2]]))
else if (identical(toupper(ylim), 'MAX')) ylim<-sapply(XY, function(x) max(x[[2]]))
else {
	ylim<-as.numeric(ylim);
	if (length(ylim[!is.na(ylim)])!=length(depth)) ylim<-sapply(XY, function(x) mean(x[[2]])+3*sd(x[[2]]));
}

if (same.ylim) ylim<-rep(max(ylim), length(depth));

WIDTH<-12;
HEIGHT<-min(18, 2*length(depth));

xlim=max(sapply(depth, length))/1000000;
ytick<-pretty(c(0, max(ylim)), n=5);
ytick<-ytick[-length(ytick)];

if (pdf) pdf('Landscape.pdf', w=WIDTH, h=HEIGHT)
else quart(w=WIDTH, h=HEIGHT);

par(mar=c(0,4,0,1), oma=c(2,0,1,0), mfrow=c(length(depth), 1));

for (i in 1:length(XY)) {
plot(XY[[i]][[1]], XY[[i]][[2]], type='h', col='darkblue', ylab=NA, xaxt='n', xaxs='i', yaxt='n', yaxs='i', xlim=c(0, xlim), ylim=c(0, ylim[i]));
axis(2, las=2, at=ytick, labels=ytick);
text(0, 0.85*ylim[i], labels=names(XY)[i], col='darkred', pos=4, cex=1.5*par()$din[2]/length(XY));
}

axis(1);

if (pdf) dev.off();
}
  
XY;
}
