Default.Annotation<-c('SYMBOL', 'CHR', 'GENENAME');

# Retrieve the annotation information of an Affy platform
SC.Retrieve.Annotation<-function(chipname, fields=c()) {
# fields    the fields of the annotation to retrieve beside the default fields as specified by Default.Annotation
l<-(paste(chipname, '.db', sep='')); 
library(l, character.only=TRUE);

fields<-c(Default.Annotation[-length(Default.Annotation)], setdiff(fields, Default.Annotation), Default.Annotation[length(Default.Annotation)]);
fields<-paste(chipname, fields, sep='');

anno<-lapply(fields, function(x) as.list(get(x)));
anno<-lapply(anno, function(x) sapply(x, function(x) if(identical(NA, x)) NA else paste(x, collapse=',')));
anno<-as.data.frame(anno);
colnames(anno)<-sub(chipname, '', fields);;
rownames(anno)<-sub('_at', '', rownames(anno));

anno;

}