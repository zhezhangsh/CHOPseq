# A wrapper of the Rsamtools:bamScanner function to read all tags on given chromosome(s) and specified SAM fields
SC.BamScannerChr<-function(bamFile, fields='DEFAULT', chr=NA, output='data.frame', by.chr=FALSE, verbose=FALSE) {

# bamFile	character, path and name of bam file
# fields	character vector or preset subset name
# chr		character vec
# output	String, the object class of returned value
#		"data.frame", default, all chromosomes in one data frame
#		"GRanges", all chromosomes in one GRanges object 
#by.chr	boolean, if true, output individual chromosomes in a list

library(Rsamtools);
library(GenomicRanges);

# define field names
if (identical(fields, 'DEFAULT')) fields<-SAM.BASIC.FIELDS
#else if (identical(fields, 'SEQ')) fields<-c('seq')
#else if (identical(fields, 'MAPQ')) fields<-c('mapq')
else if (identical(fields, 'DETAIL')) fields<-c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mapq', 'cigar', 'seq')
else if (identical(fields, 'ALL')) fields<-scanBamWhat()
#else fields<-SC.SamBasicFields();

if (length(fields)<1) fields<-SC.SamBasicFields(); # use default fields if none valid field name given
#if (output=='GRanges') fields<-union(fields, SC.SamBasicFields()); # make sure fields required by GRangs will be read in

# make sure the fields are ordered as they are in SAM format
ind<-1:length(scanBamWhat());
names(ind)<-scanBamWhat();
fields<-fields[order(ind[fields])];
######## end of define field names

	########## format output of scanBAM (list of fields) as a data frame
	format2data.frame<-function(tags) {
	tags<-lapply(tags, function(field) if (class(field)=='factor') as.vector(field)
                else if (class(field)=='DNAStringSet') as.character(field)
                else if(class(field)=='PhredQuality') as.character(field)
                else field);
	tags<-as.data.frame(tags);
	#rownames(tags)<-tags[, 'qname'];
	#tags<-tags[, names(tags)!='qname'];
	tags;
	} ###### end of function

	########### convert tags from data.frame format to a GRanges object
	format2GRanges<-function(tags, len) {
	library(GenomicRanges); 
	gr<-GRanges(seqnames=as.vector(tags[['rname']]), ranges=IRanges(start=tags[['pos']], width=tags[['qwidth']]), strand=(tags[['strand']]));
	if (length(len)==1) len<-as.numeric(len);
	seqlengths(gr)<-len;
			
	# set fields other than default as metaelement 
	#meta<-setdiff(names(tags), c("qname", "rname", "pos", "qwidth", "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", "element")); # names of meta element not include forbidden or default names	
	meta<-setdiff(names(tags), SC.SamBasicFields());
	if (length(meta)>0) {
	meta<-lapply(tags[meta], function(field) if (class(field)=='factor') as.vector(field)
                else if (class(field)=='DNAStringSet') as.character(field)
                else if(class(field)=='PhredQuality') as.character(field)
                else field);
	elementMetadata(gr)<-as.data.frame(meta);
	}	

	gr; 
	} ######### end of function

	############Assign a unique id to each alignment
	assign.id<-function(tags, chr) {
	if (class(tags)=='GRanges') names(tags)<-paste(chr, 1:length(tags), sep='_')
	else if (class(tags)=='data.frame') rownames(tags)<-paste(chr, 1:nrow(tags), sep='_')
	tags;
	}

# get chromosome names and length from BAM file header
bam.chr.len<-scanBamHeader(bamFile)[[1]][[1]];

###############################################################################################################################################
###############################################################################################################################################
if (identical(NA, chr)) {
    read.chr<-names(bam.chr.len);
}
else {
# make chromosome naming consistent, without 'chr' prefix
	# find the match of a given chromosome name to the names in BAM file
	matchName<-function(chr.name, bam.chr.names) {
		if (chr.name %in% bam.chr.names) name<-chr.name
		else name<-bam.chr.names[sub('CHR', '', toupper(bam.chr.names))==sub('CHR', '', toupper(chr.name))];
		if (length(name)==0) name<-NA;
		name;
	}
    read.chr<-sapply(chr, function(chr, bam) matchName(chr, bam), bam=names(bam.chr.len));  # the parameter passed to bamScanner()
    read.chr<-read.chr[!is.na(read.chr)];
}
#names(read.chr)<-read.chr; # used for creating GRanges object
# terminate if no valid chromosome name
if (length(read.chr)<1) {
print(as.character(chr)); 
stop("None valid chromosome names were given");
}
ranges<-RangesList(lapply(read.chr, function(chr, len) IRanges(start=1, end=len[chr]), len=bam.chr.len));
names(ranges)<-read.chr;
param<-ScanBamParam(which=ranges, what=fields); 
cat('Read in field(s):', fields, 'from BAM file\n');
##################################################################
taken<-system.time(tags<-scanBam(bamFile, param=param))[3]; # read in tags from BAM
##################################################################
names(tags)<-read.chr;

# format output
#if (identical('seq', fields) | identical('mapq', fields)) tags<-lapply(tags, function(x) x[[1]])
if (output=="GRanges" & length(setdiff(SC.SamBasicFields(), fields))==0) tags<-lapply(read.chr, function(c, tags, len) format2GRanges(tags[[c]], len[c]), tags=tags, len=bam.chr.len)
else tags<-lapply(tags, format2data.frame);
names(tags)<-read.chr;

#tags<-lapply(read.chr, function(c, tags) assign.id(tags[[c]], c), tags=tags);

if (!by.chr) {

}

tags;
}
