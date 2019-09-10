# Filter tags by chromosomes
CHOPseq.FilterByChr<-function(tags, chr, field='rname') {
# tags		data.frame, matrix, or GRanges object, must include a column or field named 'mapq'
# chr		Character vector, names of chromosomes to keep
# field		Character, column name if tags if matrix or data.frame, default is 'rname', as in SAM file

library(chipseq);
# filtering
if (class(tags)=='GRanges') tags<-tags[seqnames(tags) %in% chr] 
else if (class(tags)=='data.frame' | class(tags) =='matrix') if (field %in% colnames(tags)) tags<-tags[tags[, field] %in% chr, ] else stop(paste('\"', field, '\" is not a valid column name', sep=''))
else stop(paste('class \"', class(tags), '\" not recognized', sep=''));

tags;
}
