###############################################
ConvertCov2bedGraph<-function(cov, name, path='.', chromosome=names(cov), split.chr=FALSE, db='') {
  # cov         coverage data on all chromosomes
  # name        sample name
  # path        path to output
  # chromosome  select subset of chromsomes
  # split.chr   if split chromosomes into separate files
  # db          genome version
  
  require(IRanges);
  require(R.utils);
  
  cov<-cov[names(cov) %in% chromosome];
  
  chr<-rep(names(cov), sapply(cov, function(x) length(x@values))); 
  start<-unlist(sapply(cov, start)); 
  end<-unlist(sapply(cov, end));
  value<-unlist(sapply(cov, function(x) x@values))
  
  bed<-data.frame(chr, start, end, value, stringsAsFactors = FALSE); 

  if (split.chr) {
    file.name<-sapply(names(cov), function(chr.name) {
      header=GenerateTrackLine(type='bedGraph', name=paste(name, chr.name, sep='_'), description=paste('Depth_', name, sep=''), db=db);
      file.name<-paste('Depth_', name, '_', chr.name, '.bed', sep='');
      if(file.exists(file.name)) file.remove(file.name); 
      writeLines(header, file.name);
      write.table(bed[bed[,1]==chr.name, , drop=FALSE],  file.name, sep='\t', qu=FALSE, row=FALSE, col=FALSE, append=TRUE);
      as.vector(gzip(file.name, overwrite=TRUE));
    })
  } else {
    header=GenerateTrackLine(type='bedGraph', name=name, description=paste('Depth_', name, sep=''), db=db);
    file.name<-paste('Depth_', name, '.bed', sep='');
    if(file.exists(file.name)) file.remove(file.name); 
    writeLines(header, file.name);
    write.table(bed,  file.name, sep='\t', qu=FALSE, row=FALSE, col=FALSE, append=TRUE);
    file.name<-as.vector(gzip(file.name, overwrite=TRUE));
  }
  file.name;
}
########################

#########################################
# Generate a track line for custom track
GenerateTrackLine<-function(type='', name='', description='', visibility='', color='', itemRgb='', colorByStrand='', 
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