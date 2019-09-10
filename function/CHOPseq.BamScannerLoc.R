# A wrapper of the Rsamtools:bamScanner function to read all tags between given location on a chromosome and specified SAM fields
# Scan a bam file to read in specified fields of tags in specified region

CHOPseq.BamScannerLoc<-function(bamFile, fields='DEFAULT', chr, start=0, end=0, output='data.frame') {

# bamFile	character, path and name of bam file
# fields	character vector or preset subset name
#		DEFAULT, c('rname', 'strand', 'pos', 'qwidth', 'mapq')
#		DETAIL, c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mapq', 'cigar', 'seq')
# 		ALL, all fields
#      		custom: character vector input
# start		Integer, starting position of reading
# end		Integer, ending position of reading
# chr		character, chromosome name
# output	String, the object class of returned value, such as 'data.frame' and 'GRanges'

library(Rsamtools);
library(chipseq);

# define field names
if (fields=='DEFAULT') fields<-c('rname', 'strand', 'pos', 'qwidth', 'mapq')
else if (fields=='DETAIL') fields<-c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mapq', 'cigar', 'seq')
else if (fields=='ALL') fields<-scanBamWhat()
else fields<-intersect(fields, scanBamWhat());

if (length(fields)<1) fields<-c('rname', 'strand', 'pos', 'qwidth', 'mapq'); # use default fields if none valid field name given
if (output=='GRangs') fields<-union(fields, c('rname', 'pos', 'qwidth', 'strand')); # make sure fields required by GRangs will be read in

# make sure the fields are ordered as they are in SAM format
ind<-1:length(scanBamWhat());
names(ind)<-scanBamWhat();
fields<-fields[order(ind[fields])];
######################end of define field names

	########## format output of scanBAM (list of fields) as a data frame
	format2data.frame<-function(tags) {
	tags<-lapply(tags, function(field) if (class(field)=='factor') as.vector(field)
                else if (class(field)=='DNAStringSet') as.character(field)
                else if(class(field)=='PhredQuality') as.character(field)
                else field);
	as.data.frame(tags);
	} ###### end of function

chr<-chr[1];

# get chromosome names and length from BAM file header
bam.chr.len<-scanBamHeader(bamFile)[[1]][[1]];

# make chromosome naming consistent, without 'chr' prefix
bam.chr.names<-names(bam.chr.len);
names(bam.chr.names)<-names(bam.chr.len)<-sub('CHR', '', toupper(bam.chr.names)); 
names(chr)<-sub('CHR', '', toupper(chr));
read.chr<-bam.chr.names[sub('CHR', '', toupper(chr))]; 

if (is.na(read.chr)) {
print(chr); 
stop("Invalid choromosome name");
}

names(read.chr)<-chr;

len<-bam.chr.len[read.chr];
if (start<1|start>len) start<-1;
if (end<1|end>len) end<-len; # if start and end not specified, read in the whole chromosome
end<-pmax(start, end); # end position cannot be smaller than start

param<-ScanBamParam(which=GRanges(read.chr, IRanges(start, end)), what=fields);
tags<-scanBam(bamFile, param=param)[[1]];
tags<-format2data.frame(tags);
#rownames(tags)<-1:nrow(tags);

if ('rname' %in% fields & read.chr!=chr) tags$rname<-chr;
### end of reading

n<-nrow(tags);

# convert tags to a GRanges object
if (output=='GRanges') {
library(GenomicRanges); 
gr<-GRanges(seqnames=as.vector(tags[, 'rname']), ranges=IRanges(start=tags[, 'pos'], end=tags[, 'pos']+tags[, 'qwidth']-1), strand=as.vector(tags[, 'strand'])); 
len<-bam.chr.len[read.chr];
names(len)<-chr;
seqlengths(gr)<-len; 
meta<-setdiff(fields, c('rname', 'pos', 'qwidth', 'strand')); 
elementMetadata(gr)<-tags[, meta];
names(elementMetadata(gr))<-meta;
tags<-gr;
n<-length(tags);
}

print(paste('Retrieved', n, 'tags from', chr));
tags;
}
