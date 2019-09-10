SC.Reads2bedGraph<-function(gr, name='', by.strand=F, split.chr=F, db='') {
wd<-getwd();
if (!file.exists(name)) dir.create(name);
setwd(name);

if (!by.strand) SC.Coverage2bedGraph(coverge(gr), name, split.chr=split.chr, db=db) # output combined coverage only
else { # output coverage of each strand and their combination
if (!file.exists('Combined')) dir.create('Combined');
if (!file.exists('Forward')) dir.create('Forward');
if (!file.exists('Reverse')) dir.create('Reverse');

setwd('Combined');
SC.Coverage2bedGraph(coverage(gr), paste(name, '_Combined', sep=''), split.chr=split.chr, db=db);
setwd('../Forward');
SC.Coverage2bedGraph(coverage(gr[strand(gr)=='+']), paste(name, '_Forward', sep=''), split.chr=split.chr, db=db);
setwd('../Reverse');
SC.Coverage2bedGraph(coverage(gr[strand(gr)=='-']), paste(name, '_Reverse', sep=''), split.chr=split.chr, db=db);
}

setwd(wd);

name;
}

###############################################
SC.Coverage2bedGraph<-function(cov, name='', split.chr=F, db='') {
# cov, coverage data on all chromosomes
# name, sample name
# split.chr, if split chromosomes into separate files
# db, genome version

library(chipseq);

chr<-rep(names(cov), sapply(cov, function(x) length(x@values))); 
start<-unlist(sapply(cov, start)); 
end<-unlist(sapply(cov, end));
value<-unlist(sapply(cov, function(x) x@values))

bed<-data.frame(chr, start, end, value); 
file.name<-''; 

if (split.chr) {
write.chr<-function(chr.name, bed, name, db) {
header=SC.BedTrackLine(type='bedGraph', name=paste(name, chr.name, sep='_'), description=paste('Depth_', name, sep=''),	db=db); 
file.name<-paste('Depth_', name, '_', chr.name, '.bed', sep='');
write.table(header, file.name, sep='\t', qu=F, row=F, col=F);
write.table(bed[bed[,1]==chr.name,], file.name, sep='\t', qu=F, row=F, col=F, append=T);
SC.Gzip(file.name);
} # end of write.chr

sapply(names(cov), function(chr.name, bed, name, db) write.chr(chr.name, bed, name, db), bed=bed, name=name, db=db);
}
else {
header=SC.BedTrackLine(type='bedGraph', name=name, description=paste('Depth_', name, sep=''), db=db);
file.name<-paste('Depth_', name, '.bed', sep='');
write.table(header, file.name, sep='\t', qu=F, row=F, col=F);
write.table(bed,  file.name, sep='\t', qu=F, row=F, col=F, append=T);
SC.Gzip(file.name);
}
file.name;
}
########################



#############################################
SC.Peaks2bedGraph<-function(peaks, adjust.width=T, peak.width=1, chr.length=0, name='', desc='', db='', useScore=1, fields=0) {

peaks<-lapply(peaks, function(x) x[!is.na(x[,1]), ]);

chr<-rep(names(peaks), sapply(peaks, nrow));

peaks<-do.call(rbind, peaks);
depth<-peaks[,3]; 
center<-peaks[,1]/2+peaks[,2]/2;

if (adjust.width) {start<-floor(center-peak.width/2); end<-ceiling(center+peak.width/2);}

if (chr.length[1]!=0) {max<-chr.length[chr]; end[end>max]<-max[end>max];}

bed<-data.frame(chr, format(start, scientific=F), format(end, scientific=F))

if (fields!=0) {
id<-paste('P_', 1:length(chr), sep=''); 
bed<-cbind(bed, id, depth);
}

header<-SC.BedTrackLine(name=name, description=desc, db=db, useScore=useScore);

file<-paste('Peaks_', name, '.bed', sep='');
write.table(header, file, sep='\t', qu=F, row=F, col=F);
write.table(bed, file, sep='\t', qu=F, row=F, col=F, append=T);
SC.Gzip(file);
name;
}

#######################################################################################################
# retrieve reads from one or more samples at a specified region, convert to depth, and write to a single file
SC.Region2bedGraph<-function(gr, chr, start, end, region.name,	db='', resize=0) {
# gr                A list of named GRanges objects, each will be written out as a track
# chr, start, end   Location of the region
# region.name       Used as name of BED file
# db                Name of the genome
# resize            If greater than 0, resize GRanges objects to given number

if (resize>0) gr<-lapply(gr, function(x) resize(x, resize));

file<-paste(region.name, '.bed', sep='');

header<-paste('browser position ', chr, ':', start, '-', end, sep='');
write.table(header, file, sep='', quote=FALSE, row=FALSE, col=FALSE);

lapply(names(path), function(nm) {
x<-gr[[nm]];
x<-x[seqnames(x)==chr];
if (resize>0) x<-resize(x, resize);
cov<-coverage(x)[[chr]];
runValue(cov)[end(cov)<start]<-0;
runValue(cov)[start(cov)>end]<-0;

header=SC.BedTrackLine(type='bedGraph', name=nm, description=nm, visibility=2, db=db);
write.table(header,  file, sep='\t', qu=F, row=F, col=F, append=T);

ch<-rep(chr, nrun(cov)); 
st<-start(cov);
ed<-end(cov);
vl<-runValue(cov);
bed<-data.frame(ch, st, ed, vl); 
write.table(bed,  file, sep='\t', qu=F, row=F, col=F, append=T);

});

file;
}



#########################################
# Generate a track line for custom track
SC.BedTrackLine<-function(type='', name='', description='', visibility='', color='', itemRgb='', colorByStrand='', 
useScore='', group='', priority='', db='', offset='', url='', htmlUrl='') {

values<-c(type, name, description, visibility, color, itemRgb, colorByStrand, useScore, group, priority, db, offset, url, htmlUrl); 
names<-c('type', 'name', 'description', 'visibility', 'color', 'itemRgb', 'colorByStrand', 'useScore', 'group', 'priority', 'db', offset, 'url', 'htmlUrl');
names(values)<-names;
values[!is.na(values)&values!='']->values; 
values<-sapply(values, function(x) paste("'", x, "'", sep=''));

string<-paste(names(values), '=', values, sep='');
line<-paste(string, collapse=' ');
line<-paste('track', line, sep=' ');
line;
}
##################
