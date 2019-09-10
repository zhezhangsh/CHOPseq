# A robust wrapper of the Rsamtools:bamScanner function, allow customized reading of BAM files with a user-friendly interface
# Scan a bam file to read in specified fields of tags in specified region

CHOPseq.BamScanner<-function(bamFile, fields='DEFAULT', chr=NA, start=NA, end=NA, output='data.frame') {

# bamFile	character, path and name of bam file
# fields	character vector or preset subset name
#		DEFAULT, c('rname', 'strand', 'pos', 'qwidth', 'mapq')
#		DETAIL, c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mapq', 'cigar', 'seq')
# 		ALL, all fields
#      		custom: character vector input
# chr		character vec
# start		Integer, when all.tag is FALSE, the first base of the region to retrieve tags; the beginning of chromosome by default
# end		Integer, when all.tag is FALSE, the last base of the region to retrieve tags; the end of chromosome by default
# output	String, the object class of returned value, such as 'data.frame' and 'GRanges'

# define field names
if (fields=='DEFAULT') fields<-c('rname', 'strand', 'pos', 'qwidth', 'mapq')
else if (fields=='DETAIL') fields<-c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mapq', 'cigar', 'seq')
else if (fields=='ALL') fields<-scanBamWhat()
else fields<-intersect(fields, scanBamWhat());

if (output=='GRangs') fields<-union(fields, c('rname', 'pos', 'qwidth', 'strand')); # make sure fields required by GRangs will be read in
# make sure the fields are ordered as they are in SAM format
ind<-1:length(scanBamWhat());
names(ind)<-scanBamWhat();
fields<-fields[order(ind[fields])];

# set up parameters
if (length(chr[!is.na(chr)&chr!=''])==0) param<-ScanBamParam(what=fields) # read in all if chromosome not specified
else {
# get chromosome names and length from BAM file header
bam.chr.len<-scanBamHeader(bamFile)[[1]][[1]];
bam.chr.names<-names(bam.chr.len);

# make chromosome naming consistent, without 'chr' prefix
names(bam.chr.names)<-sub('CHR', '', toupper(bam.chr.names)); 
names(chr)<-sub('CHR', '', toupper(chr));

chr.in<-bam.chr.names[names(table(names(chr)))]; # chromosome(s) to be read in
chr.in<-chr.in[!is.na(chr.in)];

# terminate 
if (length(chr.in)<1) {
print(chr); 
stop("None valid chromosome names were given");
}


}
### end of set up parameters


# get chromosome length
#bam.trimmed.chr<-sub('CHR', '', toupper(bam.chr.names)); # make sure chromosome names are consistent (e.g. 1-22, X, Y)
#trimmed.chr<-sub('CHR', '', toupper(chr))[1]; # make sure given chromosome names are trimmed
#chr<-chr.names[which(trimmed.names==chr)];
#if (is.na(chr)) chr=chr.names[1]; # if chrosomsome name is not properly given, use the first chromosome 
#LEN<-chr.len[chr]; 

#if (is.na(start)|start>LEN) start<-0;
#if (is.na(end)|end>LEN) end<-LEN; # by default read in the whole chromosome
#end<-pmax(start, end); # end cannot be smaller than start

#param<-ScanBamParam(which=GRanges(chr, IRanges(start[1], end[1])), what=fields);
# param<-ScanBamParam(what=fields);


#tags<-scanBam(bamFile, param=param); 
#tags<-lapply(tags, function(field) if (class(field)=='factor') as.vector(field)
#               else if (class(field)=='DNAStringSet') as.character(field)
#               else if(class(field)=='PhredQuality') as.character(field)
#	       else field);
#tags<-as.data.frame(tags);
#names(tags)<-fields;

# filter reads to remove those on pseudo chromosomes
#if (all.chr&chr.num!=0) {
#if (chr.num<1|chr.num>length(chr.names)) chr.num<-length(chr.names); # get tags from all chromosomes
#chr<-chr.names[1:chr.num];
#print(chr);  
#if ('rname' %in% fields) {
#rnames<-tags[, 'rname'];
#tags<-tags[rnames %in% chr, ];
#}
#}


if (output=='GRanges') {
library(GenomicRanges); 
gr<-GRanges(seqnames=as.vector(tags[, 'rname']), ranges=IRanges(start=tags[, 'pos'], end=tags[, 'pos']+tags[, 'qwidth']-1), strand=as.vector(tags[, 'strand'])); 
seqlengths(gr)<-chr.len[names(seqlengths(gr))]; 
meta<-setdiff(fields, c('rname', 'pos', 'qwidth', 'strand')); 
elementMetadata(gr)<-tags[, meta];
names(elementMetadata(gr))<-meta;
tags<-gr;
}

#print(paste('Retrieved', nrow(tags), 'tags'));
print(chr.in);
#tags;
}
