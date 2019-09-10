# based on code of Nicola Sartori
SC.Sweatex<-function(filename,extension='Rnw',command='pdflatex',silent=FALSE,preview=FALSE) {
if (command=='latex') command='simpdftex latex --maxpfb'; 

extension<-paste('.',extension,sep='');
path=options('latexcmd')[[1]];
path=substr(path,start=1,stop=nchar(path)-5);
Sweave(paste(filename,extension,sep=''));
system(paste('./', command,' ',filename,sep=''),intern=silent);

if (preview) {
system(paste(options('pdfviewer')[[1]],' ',filename,'.pdf',sep=''))
}

}