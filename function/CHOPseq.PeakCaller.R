########################
# scan the base coverage on a chromsome and identify peaks with a minimum depth
CHOPseq.PeakCaller<-function(coverage, min.depth, ws, strand="*", merge=T, verbose=F) {
# coverage  list of Rle, each element gives the depth info. of a chromosome
# min.depth Integer, minimum depth at peak summit
# ws        Integer, window size of 2 nearest peaks
# strand    Character, indicate the strand
# merge	    Boolean, whether to merge nearby peaks with the same height as the last step
library(chipseq);

l<-length(coverage);
	
index<-lapply(coverage, function(coverage, depth) which(coverage>=depth), depth=min.depth);

# initiate peaks
peaks<-lapply(1:l, function(i, c, index) CHOPseq.InitiatePeaks(c[[i]], index[[i]]), c=coverage, index=index);
if (verbose) print(paste("Peaks initiated, N =", sum(sapply(peaks, nrow))));

# remove small peaks that have higher depth within 50 bases around
peaks<-lapply(1:l, function(i, c, p) CHOPseq.RemoveSmallPeaks(c[[i]], p[[i]]), c=coverage, p=peaks);
if (verbose) print(paste("Small peaks removed, N =", sum(sapply(peaks, nrow))));

# remove peaks shadowed by a higher peak within given window size
peaks<-lapply(1:l, function(i, c, p, ws) CHOPseq.RemoveShadowedPeaks(c[[i]], p[[i]], ws), c=coverage, p=peaks, ws=ws);
if (verbose) print(paste("Shadowed peaks removed, N =", sum(sapply(peaks, nrow))));

# if neighboring peaks have the same summit depth, remove the smaller one
peaks<-lapply(peaks, function(peaks, ws) CHOPseq.RemoveNarrowPeaks(peaks, ws), ws=ws);
if (verbose) print(paste("Narrower peaks removed, N =", sum(sapply(peaks, nrow))));

# if neighboring peaks have the same summit depth and width and close enough, merge them
peaks<-lapply(peaks, function(peaks, ws) CHOPseq.MergePeaks(peaks, ws), ws=ws);
if (verbose) print(paste("Peaks finalized, N =", sum(sapply(peaks, nrow))));


chr<-Rle(names(coverage), sapply(peaks, nrow)); 
peaks<-do.call(rbind, peaks);

if (identical(strand, '+') | identical(strand, 1)) str<-rep('+', nrow(peaks))
else if (identical(strand, '-') | identical(strand, -1)) str<-rep('-', nrow(peaks))
else str<-rep('*', nrow(peaks));

gr<-GRanges(seqnames=chr, ranges=IRanges(start=peaks[,1], end=peaks[,2]), strand=str);
seqlengths(gr)<-sapply(coverage, length);
elementMetadata(gr)<-peaks[-(1:2)]; 

names(gr)<-1:length(gr);
gr;
}

########################
# Among the indexed positions on a chromosome, find the start and end points of peaks
# Return peak information in a data frame
CHOPseq.InitiatePeaks<-function(coverage, index) {
# coverage, depth of all positions on a given chromsome
# index, positions to be checked

turn<-CHOPseq.IsTurn(coverage, index);
index<-index[!is.na(turn)]; # limit index to turning points only
turn<-turn[!is.na(turn)];
turn.before<-c(0, turn[-length(turn)]);
turn.after<-c(turn[-1], 0);

index<-index[turn==0|(turn==-1&turn.after==1)|(turn==1&turn.before==-1)];
turn<-turn[turn==0|(turn==-1&turn.after==1)|(turn==1&turn.before==-1)];

start<-index[turn==-1|turn==0]; 
end<-index[turn==1|turn==0];
depth<-as.vector(coverage[start]); 

data.frame(start=start, end=end, depth=depth); 
}
###############

#############################
# Check if indexed position on a chromosome is the turning point of depth
CHOPseq.IsTurn<-function(coverage, index) {
# coverage, depth of all positions on a given chromsome
# index, positions to be checked

depth<-as.vector(coverage[index]); 
depth.before<-rep(0, length(index));
depth.before[-1]<-as.vector(coverage[(index[-1])-1]);
depth.after<-rep(0, length(index));
depth.after[-length(index)]<-as.vector(coverage[(index[-length(index)])+1]); 

turn<-rep(NA, length(index)); # NA, not turning point; -1, turn to higher; 1, turn to lower; 0 both
turn[depth>depth.before]<--1;
turn[depth>depth.after]<-1;
turn[depth>depth.before&depth>depth.after]<-0;

turn;
}
#################


##########################
# Remove small peaks (peaks having nearby base(s) with higher depth)
CHOPseq.RemoveSmallPeaks<-function(coverage, peaks, ws=25) {
# coverage, depth at all locations on a chromosome
# peaks, a data frame of peak information, the first 3 columns are start, end and depth
# ws, window size before the start point and after the end point to search for a higher depth
	
# before<-peaks[,1]-ws;
# after<-peaks[,2]+ws;
# before[before<1]<-1;
# after[after>length(coverage)]<-length(coverage);
# max<-aggregate(coverage, start=before, end=after, FUN=max); # max depth in the region
# peaks[peaks[,3]==max, ];
	
	peaks<-peaks[order(peaks[,1]), ];
	
	depth<-peaks[,3];
	isMax<-rep(T, length(depth));
	
	for (i in 1:ws) {
		before<-peaks[,1]-i; 
		before[before<1]<-1; 
		isMax[as.vector(coverage[before])>depth]<-F;  
		
		after<-peaks[,2]+i; 
		after[after>length(coverage)]<-length(coverage); 
		isMax[as.vector(coverage[after])>depth]<-F; 
	}
	
	peaks[isMax, ];
}
##########


##############################
# Remove peaks shadowed by a nearby peak with higher depth
CHOPseq.RemoveShadowedPeaks<-function(coverage, peaks, ws) {
# coverage, depth at all locations on a chromosome
# peaks, a data frame of peak information, the first 3 columns are start, end and depth
# ws, window size of the local region to search for a higher peak

peaks<-peaks[order(peaks[,1]), ];
center<-peaks[,1]/2+peaks[,2]/2;
dist<-center[-1]-center[-length(center)]; # distance between 2 peaks
dist.before<-c(ws+1, dist); # distance to previous peak
dist.after<-c(dist, ws+1); # distance to next peak

peak0<-peaks[dist.before>ws&dist.after>ws,]; # no more work needed to these peaks
peaks<-peaks[dist.before<=ws|dist.after<=ws,];

# set non-peak regions to 0 depth
c<-coverage;
start<-c(1, peaks[,2]+1); start[start>length(c)]<-length(c); 
end<-c(peaks[,1]-1, length(c)); end[end<1]<-1;
seqselect(c, start=start, end=end)<-0;

# if the peak is the highest in the region, keep it
before<-peaks[,1]-ws;
after<-peaks[,2]+ws;
before[before<1]<-1;
after[after>length(coverage)]<-length(coverage);
max<-aggregate(coverage, start=before, end=after, FUN=max); # max depth in the region
peaks<-peaks[peaks[,3]==max, ];

peaks<-rbind(peak0, peaks);
peaks[order(peaks[,1]), ];
}
################


############################
# Remove peaks having much smaller size than its neighboring peaks
CHOPseq.RemoveNarrowPeaks<-function(peaks, ws, ratio=1/3) {
# peaks, data frame with at least 3 columns, column 1: start position; column 2: end position; column 3: depth
# ws, window size of the region to search 
# ratio, the cutoff ratio to the neighboring peaks, remove the peak if the ratio is smaller than the cutoff
	peaks<-peaks[order(peaks[,1]), ];
	
	if (nrow(peaks)>1) {
		center<-peaks[,1]/2+peaks[,2]/2;
		dist<-center[-1]-center[-length(center)];
		dist.before<-c(ws+1, dist);
		dist.after<-c(dist, ws+1);
		
		width<-peaks[,2]-peaks[,1]+1;
		width.before<-c(1, width[-length(width)]);
		width.after<-c(width[-1], 1); 
		
		temp<-peaks[(dist.before>ws|(width/width.before)>=ratio)&(dist.after>ws|(width/width.after)>=ratio), ];
		
		while(nrow(temp)!=nrow(peaks)) {
			peaks<-temp;
			
			center<-peaks[,1]/2+peaks[,2]/2;
			dist<-center[-1]-center[-length(center)];
			dist.before<-c(ws+1, dist);
			dist.after<-c(dist, ws+1);
			
			width<-peaks[,2]-peaks[,1]+1;
			width.before<-c(1, width[-length(width)]);
			width.after<-c(width[-1], 1); 
			
			temp<-peaks[(dist.before>ws|(width/width.before)>=ratio)&(dist.after>ws|(width/width.after)>=ratio), ];
		} # end of while
	} # end of if
	peaks;
}
######

###########################
# Merge peaks whose distance between center positions is smaller than given window size
CHOPseq.MergePeaks<-function(peaks, ws) {
# peaks, data frame with at least 2 columns, column 1: start position; column 2: end position
# ws, minimum distance allowed between 2 peaks
if (nrow(peaks)>1) {
peaks<-peaks[order(peaks[,1]), ];
center<-peaks[,1]/2+peaks[,2]/2;
dist<-center[-1]-center[-length(center)]; # distance between 2 peaks
dist.before<-c(ws+1, dist); # distance to previous peak
dist.after<-c(dist, ws+1); # distance to next peak

peak0<-peaks[dist.before>ws&dist.after>ws,]; # no more work needed to these peaks
peaks<-peaks[dist.before<=ws|dist.after<=ws,];

center<-peaks[,1]/2+peaks[,2]/2;
if (nrow(peaks)>1) dist<-center[-1]-center[-length(center)];

while(nrow(peaks)>1&min(dist)<=ws) {
	which<-which(dist==min(dist)); # pair of peaks having the shortest distance
	 # if there are multiple pairs of peaks having the shortest distance, randomly pick one
	if (length(which)>1) which<-sample(which, 1);

	peaks[which+1, 1]<-peaks[which, 1]; 
	peaks<-peaks[-which,];
	center<-peaks[,1]/2+peaks[,2]/2;
	if (nrow(peaks)>1) dist<-center[-1]-center[-length(center)]; 

} # end of while

peaks<-rbind(peak0, peaks);
peaks<-peaks[order(peaks[,1]),];		
} # end of if
peaks;
}
##########
