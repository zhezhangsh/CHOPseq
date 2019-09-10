SC.DAVIDannots<-list(
	All=data.frame(
	Category=c('Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Pathways', 'Pathways', 'Pathways', 'Pathways', 'Pathways', 'Pathways', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'General Annotations', 'Functional Categories', 'Functional Categories', 'Functional Categories', 'Functional Categories', 'Functional Categories', 'Functional Categories', 'Functional Categories', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Protein-Protein Interaction', 'Literature', 'Literature', 'Literature', 'Disease', 'Disease'), 
	Term=c('GOTERM_BP_1', 'GOTERM_BP_2', 'GOTERM_BP_3', 'GOTERM_BP_4', 'GOTERM_BP_5', 'GOTERM_BP_ALL', 'GOTERM_BP_FAT', 'GOTERM_CC_1', 'GOTERM_CC_2', 'GOTERM_CC_3', 'GOTERM_CC_4', 'GOTERM_CC_5', 'GOTERM_CC_ALL', 'GOTERM_CC_FAT', 'GOTERM_MF_1', 'GOTERM_MF_2', 'GOTERM_MF_3', 'GOTERM_MF_4', 'GOTERM_MF_5', 'GOTERM_MF_ALL', 'GOTERM_MF_FAT', 'BLOCKS_ID', 'COG', 'INTERPRO', 'PDB_ID', 'PFAM', 'PIR_ALN', 'PIR_HOMOLOGY_DOMAIN', 'PIR_SUPERFAMILY', 'PRINTS', 'PRODOM', 'PROSITE', 'SCOP_ID', 'SMART', 'TIGRFAMS', 'BBID', 'BIOCARTA', 'EC_NUMBER', 'KEGG_COMPOUND', 'KEGG_PATHWAY', 'KEGG_REACTION', 'ALIAS_GENE_SYMBOL', 'CHROMOSOME', 'CYTOBAND', 'GENE', 'GENE_SYMBOL', 'HOMOLOGOUS_GENE', 'LL_SUMMARY', 'OMIM_ID', 'PIR_SUMMARY', 'PROTEIN_MW', 'REFSEQ_PRODUCT', 'SEQUENCE_LENGTH', 'SP_COMMENT', 'CGAP_EST_QUARTILE', 'CGAP_EST_RANK', 'COG_ONTOLOGY', 'PIR_SEQ_FEATURE', 'SP_COMMENT_TYPE', 'SP_PIR_KEYWORDS', 'UP_SEQ_FEATURE', 'BIND', 'DIP', 'HIV_INTERACTION_CATEGORY', 'HIV_INTERACTION', 'MINT', 'NCICB_CAPATHWAY', 'TRANSFAC_ID', 'GENERIF_SUMMARY', 'HIV_INTERACTION_PUBMED_ID', 'PUBMED_ID', 'GENETIC_ASSOCIATION_DB_DISEASE', 'OMIM_DISEASE')),
	
	Favorites=data.frame(
	Category=c('Gene_Ontology', 'Gene_Ontology', 'Gene_Ontology', 'Protein_Domains', 'Protein_Domains', 'Protein_Domains', 'Protein-Protein Interaction', 'Pathways', 'Pathways'),
	Term=c('GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', 'INTERPRO', 'PIR_SUPERFAMILY', 'SMART', 'TRANSFAC_ID', 'BIOCARTA', 'KEGG_PATHWAY'))
);

SC.DAVIDtools<-data.frame(
	DAVID.function=c('Gene Functional Classification', 'Funtional Annotation Clustering', 'Functional Annotation Summary', 'Functional Annotation Chart', 'Functional Annotation Table', 'Show Gene List Names in Batch', 'Gene Report', 'Gene Full Report'), 
	Valid.value=c('gene2gene ', 'term2term', 'summary', 'chartReport', 'annotationReport', 'list', 'geneReport', 'geneReportFull')
);

SC.DAVIDtypes<-c('AFFYMETRIX_3PRIME_IVT_ID', 'AFFYMETRIX_EXON_GENE_ID', 'AFFYMETRIX_SNP_ID', 'AGILENT_CHIP_ID', 'AGILENT_ID', 'AGILENT_OLIGO_ID', 'ENSEMBL_GENE_ID', 'ENSEMBL_TRANSCRIPT_ID', 'ENTREZ_GENE_ID', 'FLYBASE_GENE_ID', 'FLYBASE_TRANSCRIPT_ID', 'GENBANK_ACCESSION', 'GENPEPT_ACCESSION', 'GENOMIC_GI_ACCESSION', 'PROTEIN_GI_ACCESSION', 'ILLUMINA_ID', 'IPI_ID', 'MGI_ID', 'GENE_SYMBOL', 'PFAM_ID', 'PIR_ACCESSION', 'PIR_ID', 'PIR_NREF_ID', 'REFSEQ_GENOMIC', 'REFSEQ_MRNA', 'REFSEQ_PROTEIN', 'REFSEQ_RNA', 'RGD_ID', 'SGD_ID', 'TAIR_ID', 'UCSC_GENE_ID', 'UNIGENE', 'UNIPROT_ACCESSION', 'UNIPROT_ID', 'UNIREF100_ID', 'WORMBASE_GENE_ID', 'WORMPEP_ID', 'ZFIN_ID');


# Perform a DAVID "Functional Annotation Chart" analysis 
SC.DAVIDchart<-function(ids, type=SC.DAVIDtypes[9], annot=as.vector(SC.DAVIDannots$Favorites[,2]), use.symbol=TRUE, gene.symbol=c()) {
# ids		Gene IDs for query
# type		Gene ID type, default is Entrez gene ID
# annot		Types of gene class terms, such as GO_BP
# use.symbol	Whether to use gene symbols in the result table. If true and gene.symbol not properly set, query DAVID to map gene symbols
# gene.symbol   If named by ids and use.symbol is TRUE, use gene.symbol to id mapping information in this field in DAVID result 

library(DAVIDQuery);

#if (type=='REFSEQ_MRNA') ids<-ids[grep('NM_', ids)];
mapping<-NULL;
try (mapping<-DAVIDQuery(ids, type=type, tool='geneIdConversion', annot='DAVID', details=FALSE), silent=TRUE);

if (is.null(mapping)) {
results<-as.data.frame(matrix(rep('', 8), nr=1, nc=8));
names(results)<-c('Category', 'Term', 'Count', '%', 'Fold Enrichment', 'PValue', 'FDR', 'Genes');
} 
else {
mapped<-split(as.vector(mapping[,1]), mapping[,2])

options(warn=-1);
cat('\n', '##############', '\n', 'Please ignore the following error message(s)', '\n');
results<-lapply(annot, function(a, ids) try(DAVIDQuery(ids=ids, type='DAVID', tool='chartReport', annot=a, details=FALSE), silent=TRUE), ids=names(mapped)); 
options(warn=1);
cat('##############', '\n\n');

results<-results[sapply(results, class)!='try-error']; # if too few IDs are recognized, no results will be returned.

if (length(results)==0) {
results<-as.data.frame(matrix(rep('', 8), nr=1, nc=8));
names(results)<-c('Category', 'Term', 'Count', '%', 'Fold Enrichment', 'PValue', 'FDR', 'Genes');
}
else { 
results<-lapply(results, function(x) x[-1,]);
results<-do.call('rbind', results);
results<-results[, c(1:4, 10, 5, 13, 6)];
for (i in 3:7) results[[i]]<-as.numeric(as.vector(results[[i]]));
results<-results[order(results[,6]), ];

# map back
x<-strsplit(as.vector(results[,ncol(results)]), ', ');  
x<-lapply(x, function(x, mapped) unlist(mapped[x]), mapped=mapped);
x<-lapply(x, function(x) names(table(x))); 
x<-sapply(x, function(x) paste(x, collapse=', ')); 
results[, ncol(results)]<-x; 

# use gene symbol in the result reports
if (use.symbol) { # use gene symbol in report
if (identical(ids, names(gene.symbol[ids]))) { # use gene symbol in report
x<-strsplit(as.vector(results[,ncol(results)]), ', ');  
x<-lapply(x, function(x, smb) as.vector(smb[x]), smb=gene.symbol);
x<-sapply(x, function(x) paste(x, collapse=', ')); 
results[, ncol(results)]<-x;   
} # end of 'if (identi ...'
else {
mapping<-NULL;
try (mapping<-DAVIDQuery(ids, type=type, tool='geneIdConversion', annot='OFFICIAL_GENE_SYMBOL', details=FALSE), silent=TRUE);
if (!is.null(mapping)) { 
    mapped<-split(mapping[,2], mapping[,1]);
    x<-strsplit(as.vector(results[,ncol(results)]), ', ');  
    x<-lapply(x, function(x, mapped) unlist(mapped[x]), mapped=mapped); 
    x<-lapply(x, function(x) names(table(x))); 
    x<-sapply(x, function(x) paste(x, collapse=', ')); 
    results[, ncol(results)]<-x; 
} # end of 'if (!is.null(mapping))'  
} # end of else
} # end of use.symbol
} # end of length(results)>0
} # end of !is.null(mapping)


results;
}

