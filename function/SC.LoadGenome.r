# Load BSgenome library of a given genome
SC.LoadGenome<-function(genome.name) {
# genome.name   Example: 'hg19', 'human', 'Home sapien'

g<-tolower(genome.name); 
  
if (g %in% c('hg19', 'human', 'homo sapiens', 'hs', 'h.s.', 'hsapiens', 'grch37')) {
  library(BSgenome.Hsapiens.UCSC.hg19);
  x<-'Hsapiens';
}
else if (g %in% c('hg18', 'ncbi36')) {
  library(BSgenome.Hsapiens.UCSC.hg18);
  x<-'Hsapiens';
}
else if (g %in% c('mm9', 'mouse', 'mouse musculus', 'mmusculus', 'ncbi37')) {
  library(BSgenome.Mmusculus.UCSC.mm9);
  x<-'Mmusculus';
}
else if (g %in% c('mm8')) {
  library(BSgenome.Mmusculus.UCSC.mm8);
  x<-'Mmusculus';
}
else stop('Error: Unknown genome name', genome.name)

x;
}
