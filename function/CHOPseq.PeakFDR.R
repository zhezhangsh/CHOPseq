# Calculate the FDR of peaks based on their number of reads in 2 compared samples
CHOPseq.PeakFDR<-function(read1, read0, adjust=T, cutoff.count=c(0, Inf), cutoff.fc=c(0, Inf), select=0) {
# read1: a matrix or data frame with 2 columns; each row represents a peak from ChIP sample
# read0: a matrix or data frame with 2 columns; each row represents a peak from control sample
 # the first column is the read count of peaks in the sample from which the peaks were identified, 
 # the second column is the read count of the peak regions in the other sample
# adjust: whether adjust FDR by total number of peaks from both samples
# cutoff.count & cutoff.fc: the cutoff values of read count and fold change, used to select a subset of peaks
 # when the number of peaks is large, the running time will be long, raise cutoffs to reduce running time and improve specificity
# select: select given number of control peaks if there are more peaks; raise this value may increase running time significantly,
 # no selection if equal to 0
# return a numeric vector the same length as the number of ChIP peaks
data.frame(read1)->read1;
data.frame(read0)->read0; 

read1[read1[,2]<1,2]<-1; # round to 1 if there is no read in the peak region
read0[read0[,2]<1,2]<-1; 

read1[,2]<-read1[,1]/read1[,2]; # fold enrichment of reads at all peak regions
read0[,2]<-read0[,1]/read0[,2]; 

fdr<-rep(NA, nrow(read1)); 
names(fdr)<-rownames(read1)<-1:nrow(read1); # keep the indexe of peaks after the filtering of peaks in the next step
read1<-read1[read1[,1]>=cutoff.count[1]&read1[,1]<=cutoff.count[2]&read1[,2]>=cutoff.fc[1]&read1[,2]<=cutoff.fc[2], ];
read0<-read0[read0[,1]>=cutoff.count[1]&read0[,1]<=cutoff.count[2]&read0[,2]>=cutoff.fc[1]&read0[,2]<=cutoff.fc[2], ];

read<-read1; # compare ChIP peaks to this object
if (select>0&nrow(read)>select) read<-read[sample(1:nrow(read), select),];  
if (select>0&nrow(read0)>select) read0<-read0[sample(1:nrow(read0), select), ];

m<-sapply(1:nrow(read1), function(i, read1, read) nrow(read[read[, 1]>=read1[i, 1]&read[, 2]>=read1[i, 2],]), read1=read1, read=read);
n<-sapply(1:nrow(read1), function(i,read1,read0) nrow(read0[read0[,1]>=read1[i,1]&read0[,2]>=read1[i,2],]), read1=read1, read0=read0);

fdr[rownames(read1)]<-n/m;
if (adjust) fdr<-fdr/(nrow(read0)/nrow(read1)); # adjust FDR based on total number of qualified peaks from the 2 samples
fdr[is.na(fdr)]<-1;
fdr[fdr>1]<-1; # FDR cannot be greater than 1

as.vector(fdr);
}
