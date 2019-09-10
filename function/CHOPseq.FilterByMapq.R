# Filter tags by their mapping quality given as the 'mapq' value in SAM files
CHOPseq.FilterByMapq<-function(tags, scores, remove=F) {
# tags		data.frame, matrix, or GRanges object, must include a column or field named 'mapq'
# scores	Integer, 'mapq' scores will allow tags to be accepted
# remove	Boolean, whether to remove the mapq information after the filtering
library(chipseq);

# filtering
if (class(tags)=='GRanges') { # if tags are in "GRanges" format, the mapq must be excluded in elementMetadata as a field names 'mapq'
meta<-elementMetadata(tags);
if ('mapq' %in% colnames(meta)) { # if mapq information exists
mapq<-meta[, 'mapq'];
if (remove) elementMetadata(tags)<-meta[, colnames(meta)[colnames(meta)!='mapq']];
tags<-tags[mapq %in% scores]; 
}
} # end of if (class(tags)=='GRanges')
# if tags are in a matrix or data.frame, the mapq must be in a field named 'mapq'
else {
if ('mapq' %in% colnames(tags)) { # if mapq information exists
mapq<-tags[, 'mapq'];
if (remove) tags<-tags[, colnames(tags)[colnames(tags)!='mapq']]; 
tags<-tags[mapq %in% scores,]
}
} # end of if/else

tags;
}
