# load the reference genome library according to genome name
SC.LoadReference<-function(genome.name) {
  if (toupper(genome.name) == 'HG19') {
    library(BSgenome.Hsapiens.UCSC.hg19);
    'Hsapiens';
  } else if (toupper(genome.name) == 'MM9') {
    library(BSgenome.Mmusculus.UCSC.mm9);
    'Mmusculus';
  } else if (toupper(genome.name) == 'MM8') {
    library(BSgenome.Mmusculus.UCSC.mm8);
    'Mmusculus';
  } else {NA} 
}