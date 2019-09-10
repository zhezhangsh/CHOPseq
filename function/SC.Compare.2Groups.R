# Perform a default analysis of 2-group comparison
SC.Compare.2Groups<-function(
	d, 
	anno, 
	ind.g1=1:round(ncol(d)/2), 
	ind.g2=(max(ind.g1)+1):ncol(d), 
	groupnames=c('A', 'B'), 
	paired=FALSE, 
	DAVID=TRUE,
	GSEA=TRUE,
	volcano=TRUE, 
	pca=TRUE,
	id.type="ENTREZ_GENE_ID", 
	GSEA.Set=list(),
	cutoff.l2r=0, 
	cutoff.p=0.05, 
	num.top=250, 
	write.out=TRUE
) {

anno<-anno[rownames(d), ];

filename<-paste(groupnames[1], groupnames[2], sep=' vs. ');
wd<-getwd();

cat('Putting results to ', file.path(wd, filename), '\n');
if (file.exists(filename)) { # already exists as a file or directory
	if (!file.info(filename)$isdir) {
		file.rename(filename, paste(filename, 'file', sep='_')); # if file with the same name exist, rename the file
		dir.create(filename);
	}
}
else dir.create(filename);

# run SAM 
cat('Running SAM ...', '\n');
sam<-CHOPbic.SAM.2Groups(d, ind.g1, ind.g2, groupnames, paired, 100);
out<-cbind(anno[rownames(sam),], sam);

up<-out[out[, 'p_SAM']<=cutoff.p&out[, ncol(anno)+3]>=cutoff.l2r, ];
dn<-out[out[, 'p_SAM']<=cutoff.p&out[, ncol(anno)+3]<=(-1*cutoff.l2r), ];
up<-up[order(up[, ncol(anno)+3], decreasing=TRUE),][1:min(num.top, nrow(up)), ];
dn<-dn[order(dn[, ncol(anno)+3]), ][1:min(num.top, nrow(dn)), ];

# parameters for formatting
f1<-list(id.name='Gene_ID', format=data.frame(ind=ncol(anno)+(1:6), fmt=c(rep('0.0000', 4), '0.0E+0', '0.00%'))); # format: group comparison
f2<-list(row.names=FALSE, format=data.frame(ind=3:7, fmt=c('0', '0.00', '0.00', '0.0E+0', '0.0E+0'))); # format: DAVID output
f3<-list(row.names=FALSE, format=data.frame(ind=4:13, fmt=rep('0.000', 10))); # format: GSEA output

# output of group comparison
settings<-list(f1, f1, f1);
values<-list(out, up, dn);
names(values)<-c('Complete List', paste('Top', nrow(up), ', higher in ', groupnames[2], sep=''), paste('Top', nrow(dn), ', lower in ', groupnames[2], sep=''));

# run DAVID analysis
if (DAVID) {
cat('Running DAVID ...', '\n');
i<-grep('symbol', colnames(anno), ignore.case=TRUE);
if (is.na(i[1])) symb<-c() else {symb<-as.vector(anno[,i[1]]); names(symb)<-rownames(anno);}
dv.up<-CHOPbic.DAVID.chart(rownames(up), type=id.type, gene.symbol=symb);
dv.dn<-CHOPbic.DAVID.chart(rownames(dn), type=id.type, gene.symbol=symb);

settings<-append(settings, list(f2, f2));
values[[length(values)+1]]<-dv.up;
names(values)[length(values)]<-paste('DAVID, higher in ', groupnames[2], sep='');

values[[length(values)+1]]<-dv.dn;
names(values)[length(values)]<-paste('DAVID, lower in ', groupnames[2], sep='');
}

if (GSEA) {
cat('Running GSEA ...', '\n');
wd<-getwd();
setwd(paste(wd, filename, sep='/'));
i<-grep('symbol', colnames(anno), ignore.case=TRUE);
if (is.na(i[1])) symb<-rownames(d) else {symb<-as.vector(anno[,i[1]]); names(symb)<-rownames(anno);}
gsea<-list();
nm<-names(GSEA.Set);
if (is.null(nm)) nm<-paste('Collection', 1:length(GSEA.Set), sep='_')
else for (i in 1:length(nm)) if (is.na(nm[i])) nm[i]<-paste('Collection', i, sep='_')
for (i in 1:length(GSEA.Set)) {
gsea.folder<-paste('GSEA-', nm[i], sep=''); 
if (!file.exists(gsea.folder)) dir.create(gsea.folder);
setwd(gsea.folder);
gsea[[i]]<-CHOPbic.GSEA(d, symb, ind.g1, ind.g2, groupnames, gs.size.cutoff = c(10, 500), gmt=GSEA.Set[[i]]);
names(gsea[[i]])[4]<-paste('ES, ', groupnames[2], '-', groupnames[1], sep='');
names(gsea[[i]])[5]<-paste('NES, ', groupnames[2], '-', groupnames[1], sep='');

if (file.exists("All_Plots.pdf")) file.copy("All_Plots.pdf", "../All_Plots.pdf");
setwd('..');
if (file.exists("All_Plots.pdf")) file.rename("All_Plots.pdf", paste(gsea.folder, "-Plots.pdf", sep=''));
}
for (i in 1:length(nm)) settings[[length(settings)+1]]<-f3;
values[(length(values)+1):(length(values)+length(nm))]<-gsea;
names(values)[(length(values)-length(nm)+1):length(values)]<-paste('GSEA,', nm);

setwd(wd);
}

# write to Excel file
if (write.out) {
cat('Writing Excel ...', '\n');
setwd(paste(wd, filename, sep='/'));
CHOPbic.Excel(values, fileName=filename, settings=settings);
setwd(wd);
}

# Plotting
if (volcano|pca) cat('Plotting ...', '\n');
if (volcano) {
    setwd(paste(wd, filename, sep='/'));
	CHOPbic.Volcano(sam[,3], sam[,5], groupnames, screen=FALSE, pdf=TRUE, tiff=TRUE);
	setwd(wd);
    cat('plotted Volcano plot\n');
}
if (pca) {
	i<-grep('chr', colnames(anno), ignore.case=TRUE);
	if (length(i)>0) {
		chr<-as.vector(anno[, i[1]]);
		xy<-union(grep('x', chr, ignore.case=TRUE), grep('y', chr, ignore.case=TRUE));
		e<-d[setdiff(1:nrow(d), xy),];
	}
	else e<-d;
	
	setwd(paste(wd, filename, sep='/'));
	SC.PCA.2Groups(prcomp(t(e[, c(ind.g1, ind.g2)])), 1:length(ind.g1), length(ind.g1)+(1:length(ind.g2)), TRUE, paste('PCA,', filename), groupnames);
	setwd(wd);
    cat('plotted PCA plot\n');
}

# Create and write out log
setwd(paste(wd, filename, sep='/'));
log<-list(data=d, anno=anno, ind.g1=ind.g1, ind.g2=ind.g2, groupnames=groupnames, paired=paired, 
    DAVID=DAVID, GSEA=GSEA, volcano=volcano, pca=pca, id.type=id.type, GSEA.Set=GSEA.Set, 
    cutoff.l2r=cutoff.l2r, cutoff.p=cutoff.p, num.top=num.top, write.out=write.out, folder=getwd());
save(log, file='log.RData');

setwd(wd);

values;
}