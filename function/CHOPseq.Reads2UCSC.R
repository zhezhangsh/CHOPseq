CHOPseq.Reads2bedGraph<-function(gr, name='', by.strand=F, split.chr=F, db='') {
wd<-getwd();
if (!file.exists(name)) dir.create(name);
setwd(name);

if (!by.strand) CHOPseq.Coverage2bedGraph(coverge(gr), name, split.chr=split.chr, db=db) # output combined coverage only
else { # output coverage of each strand and their combination
if (!file.exists('Combined')) dir.create('Combined');
if (!file.exists('Forward')) dir.create('Forward');
if (!file.exists('Reverse')) dir.create('Reverse');

setwd('Combined');
CHOPseq.Coverage2bedGraph(coverage(gr), paste(name, '_Combined', sep=''), split.chr=split.chr, db=db);
setwd('../Forward');
CHOPseq.Coverage2bedGraph(coverage(gr[strand(gr)=='+']), paste(name, '_Forward', sep=''), split.chr=split.chr, db=db);
setwd('../Reverse');
CHOPseq.Coverage2bedGraph(coverage(gr[strand(gr)=='-']), paste(name, '_Reverse', sep=''), split.chr=split.chr, db=db);
}

setwd(wd);
}

###############################################
CHOPseq.Coverage2bedGraph<-function(cov, name='', split.chr=F, db='') {
# cov, coverage data on all chromosomes
# name, sample name
# split.chr, if split chromosomes into separate files
# db, genome version

library(chipseq);

    #########################################
    # Generate a track line for custom track
    trackLine<-function(type='', name='', description='', visibility='', color='', itemRgb='', colorByStrand='', 
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

    # gzip a file
    gzip<-function(file.name, remove=T) {
    # file.name; the name of the file to be compressed
    # remove; whether keep the original file
    if (remove) system(paste('gzip -f ', file.name))
    else system(paste('gzip -c -f ', file.name, ' > ', file.name, '.gz', sep='' ));
    }

chr<-rep(names(cov), sapply(cov, function(x) length(x@values))); 
start<-unlist(sapply(cov, start)); 
end<-unlist(sapply(cov, end));
value<-unlist(sapply(cov, function(x) x@values))

bed<-data.frame(chr, start, end, value); 
file.name<-''; 

if (split.chr) {
write.chr<-function(chr.name, bed, name, db) {
header=trackLine(type='bedGraph', name=paste(name, chr.name, sep='_'), description=paste('Depth_', name, sep=''),	db=db); 
file.name<-paste('Depth_', name, '_', chr.name, '.bed', sep='');
write.table(header, file.name, sep='\t', qu=F, row=F, col=F);
write.table(bed[bed[,1]==chr.name,], file.name, sep='\t', qu=F, row=F, col=F, append=T);
gzip(file.name);
} # end of write.chr

sapply(names(cov), function(chr.name, bed, name, db) write.chr(chr.name, bed, name, db), bed=bed, name=name, db=db);
}
else {
header=trackLine(type='bedGraph', name=name, description=paste('Depth_', name, sep=''), db=db);
file.name<-paste('Depth_', name, '.bed', sep='');
write.table(header, file.name, sep='\t', qu=F, row=F, col=F);
write.table(bed,  file.name, sep='\t', qu=F, row=F, col=F, append=T);
gzip(file.name);
}
name;
}
########################



#############################################
CHOPseq.Peaks2bed<-function(peaks, adjust.width=T, peak.width=1, chr.length=0, name='', desc='', db='', useScore=1, fields=0) {


    #########################################
    # Generate a track line for custom track
    trackLine<-function(type='', name='', description='', visibility='', color='', itemRgb='', colorByStrand='', 
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

    # gzip a file
    gzip<-function(file.name, remove=T) {
    # file.name; the name of the file to be compressed
    # remove; whether keep the original file
    if (remove) system(paste('gzip -f ', file.name))
    else system(paste('gzip -c -f ', file.name, ' > ', file.name, '.gz', sep='' ));
    }

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

header<-trackLine(name=name, description=desc, db=db, useScore=useScore);

file<-paste('Peaks_', name, '.bed', sep='');
write.table(header, file, sep='\t', qu=F, row=F, col=F);
write.table(bed, file, sep='\t', qu=F, row=F, col=F, append=T);
gzip(file);
name;
}

