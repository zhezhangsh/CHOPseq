# process tags: filter, draw statistics, and output UCSC files
CHOPseq.ProcessTag<-function(path, name, mapq, resize, db='', UCSC=F, verbose=F) {
# path		Character, path and name of file saves GRanges object
# name		Character, sample name
# mapq		Integer vector, mapq scores not to be filtered
# resize	Integer, resize the tag length to
# db		Character, genome build symbol
# UCSC		Boolean, whether write coverage data for UCSC visualization

# create a folder under current working directory 
# name it  with sample name, and save all outputs there
library(chipseq);

if (!file.exists(name)) dir.create(name);
folder<-paste(getwd(), '/', name, sep='');

load(path)->temp;
gr<-get(temp);
if (verbose) print('Tags loaded');

# filtering by mapq, and save data after filtering to disk
gr.total<-gr;
keep<-elementMetadata(gr)$mapq %in% mapq;
gr.unused<-gr[!keep];
gr<-gr[keep];
save(gr.unused, file=paste(folder, '/GR_', name, '_unused.RData', sep=''));
save(gr, file=paste(folder, '/GR_', name, '_used.RData', sep=''));
if (verbose) print('Tags filtered');

gr<-resize(gr, width=resize);

# convert to coverage data and save to disk
cov<-coverage(gr);
cov.fw<-coverage(gr[strand(gr)=='+']);
cov.rv<-coverage(gr[strand(gr)=='-']);
save(cov, file=paste(folder, '/Coverage_', name, '.RData', sep=''));
save(cov.fw, file=paste(folder, '/Coverage_Forward_', name, '.RData', sep=''));
save(cov.rv, file=paste(folder, '/Coverage_Reverse_', name, '.RData', sep=''));
if (verbose) print('Converted to coverage data');

######################################## convert to UCSC bedGraph file
if (UCSC) {
wd<-getwd();
if (!file.exists('UCSC')) dir.create('UCSC');
setwd('UCSC');
if (!file.exists(name)) dir.create(name);
setwd(name);
if (!file.exists('Combined')) dir.create('Combined');
if (!file.exists('Forward')) dir.create('Forward');
if (!file.exists('Reverse')) dir.create('Reverse');

if (verbose) print('Writing UCSC bedGraph files...');
setwd('Combined');
CHOPseq.Coverage2bedGraph(cov, paste(name, '_Combined', sep=''), split.chr=T, db=db);
if (verbose) print('Combined strands');
setwd('../Forward');
CHOPseq.Coverage2bedGraph(cov.fw, paste(name, '_Forward', sep=''), split.chr=T, db=db);
if (verbose) print('Forward strand');
setwd('../Reverse');
CHOPseq.Coverage2bedGraph(cov.fw, paste(name, '_Reverse', sep=''), split.chr=T, db=db);
if (verbose) print('Reverse strand');
setwd(wd);
} # end of write UCSC

# get stats by chromosomes
len<-seqlengths(gr);
chr.names<-names(len);
eff.size<-CHOPseq.EffectiveSize(cov);
tag.used<-table(seqnames(gr));
tag.total<-table(seqnames(gr.total));
fraglen<-estimate.mean.fraglen(gr.total);
depth<-mean(cov);

chr.names<-c(chr.names, 'Overall');
len<-c(len, sum(as.numeric(len)));
tag.total<-c(tag.total, sum(tag.total));
tag.used<-c(tag.used, sum(tag.used));
fraglen<-c(fraglen, weighted.mean(fraglen, w=eff.size));
depth<-c(depth, weighted.mean(depth, w=eff.size));
eff.size<-c(eff.size, sum(as.numeric(eff.size)));
chr.stats<-data.frame(Length=len, Effective.Length=eff.size, Total.Tag=tag.total, 
Used.Tag=tag.used, Fragment.Length=fraglen, Mean.Depth=depth);
rownames(chr.stats)<-chr.names;
save(chr.stats, file=paste(folder, '/Chr_Stats_', name, '.RData', sep=''));

print(paste('Finished processing sample', name));

chr.stats;
}
