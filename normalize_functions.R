

run_deseq2 <- function(file_path, save_name, group1_terms, group1_name, group2_terms, group2_name){
library(DESeq2)
raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
counts <- raw.data[ , -c(1,ncol(raw.data)) ]
rownames(counts) <- row.names(raw.data)
group <- c()
for (x in colnames( counts )){if(grepl(group1_terms,x)){ group <- append(group, group1_name)} else if(grepl(group2_terms,x)){ group <- append(group,group2_name)} }
print(colnames(counts))
print(group)
samples <- data.frame(row.names=colnames(counts), condition=as.factor(c(group)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
matrix <- counts(dds,normalized=TRUE)
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
write.table(matrix, file=paste("DESeq_", save_name,"matrix_norm.txt", sep = "_"), sep="\t", col.names=NA)
saveRDS(dds, paste(save_name,".rds"))
plotDispEsts(dds)
return(dds)}

make_cpm <- function(path_to_file,file_base_name){
  library(edgeR)
	raw.data <- read.csv(file = file_path, header=TRUE, sep="\t", row.names=1)
	counts <- raw.data[ , -c(1,ncol(raw.data)) ]
	rownames( counts ) <- row.names(raw.data)
	cds <- DGEList( counts , group = colnames(counts) )
	cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
	cds <- calcNormFactors(cds, method="TMM")
	cps <- cpm(cds, normalized.lib.sizes=TRUE)
	write.table(cps, file=paste(file_base_name,"normalized_edgR_cpm.txt",sep='_'), sep="\t", col.names=NA)
}
