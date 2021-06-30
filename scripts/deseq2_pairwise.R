library(DESeq2)
library(rhdf5)
library("tximport")
library(GenomicFeatures)

library(mygene)
library(pheatmap)
library(genefilter)
library(ggplot2)


abundance_files<-snakemake@input[["abundance_h5"]]
dge_names<-snakemake@params[["dge_names"]]
print(dge_names)

dge_names=gsub(" ", "", dge_names, fixed = TRUE)
print(dge_names)
gtf<-snakemake@params[["gtf_AS_isoforms"]]
sample<-matrix(dge_names, ncol=1)
colnames(sample)<-"sample"
rownames(sample)=seq(1, length(dge_names))
samples<-as.data.frame(sample)
print(samples)
files=file.path(abundance_files)
names(files)=samples$sample

txdb<-makeTxDbFromGFF(gtf)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

print("revised tx2gene")
head(tx2gene)
imported_ab<-tximport(files, type="kallisto", tx2gene = tx2gene)

sampleTable<-data.frame(condition=factor(dge_names))
print(sampleTable)

rownames(sampleTable)<-colnames(imported_ab$counts)






dds<-DESeqDataSetFromTximport(imported_ab, sampleTable, ~condition)
vst_dds <- vst(dds, blind=TRUE, fitType = "local")
vst=assay(vst_dds)

#pca
print(snakemake@output[["pca"]])
par(oma=c(2,2,2,2))
png(snakemake@output[["pca"]], width = 25, height = 15, units = "cm",res=300, pointsize =12)
pca_data=plotPCA(vst_dds, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color=condition)) +
  geom_point() +  
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(legend.title = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(color = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  coord_fixed()
dev.off()

#replicates
print(length(sampleTable$condition)!=unique(sampleTable$condition))
print( (length(sampleTable$condition)) ==(unique(sampleTable$condition)))
if(any((length(sampleTable$condition))==(unique(sampleTable$condition)) ) )
{

dds <- DESeq(dds)

contrast_combinations=combn(levels(sampleTable$condition), 2, FUN=list)
result_list=list()

for(i in contrast_combinations)
{
  #print(i)
  result_list[[paste0(i[1],"_vs_", i[2])]]<-results(dds, contrast = c("condition", i[1], i[2]))
  # print(nrow(result_list[[paste0(i[1],"_vs_", i[2])]]))
  genenames_split=sapply(rownames(result_list[[paste0(i[1],"_vs_", i[2])]]), strsplit, "\\.")#[[1]][1])
  
  query_list=lapply(genenames_split, `[[`, 1)
  gene_symbols=getGenes(query_list, field="symbol")
  gene_symbols$symbol <- ifelse(is.na(gene_symbols$symbol), gene_symbols$query, gene_symbols$symbol)

  result_list[[paste0(i[1],"_vs_", i[2])]]$Gene<-gene_symbols$symbol
  
}

for(i in names(result_list))
{
  write.table(as.data.frame(result_list[[i]]), file=file.path(cwd, paste0(i,".diffgenes",".tab")), sep="\t" , quote = FALSE, row.names = FALSE)
  
}

for(i in names(result_list))
{
  
  sig_genes=append(sig_genes, rownames(subset(result_list[[i]], padj<=padj_cutoff)))
  #plot_name=(paste0(i,"_"))
  # plot_name=append(substitute(i,"_"))
}

sig_genes=unique(sig_genes)


vst_sig=vst[unlist(sig_genes), ]

#str(vst)
z_score <- function(x){(x - mean(x))/sd(x)}
as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(sampleTable$condition))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_sig)
vst_sig_z_scored <- t(apply(vst_sig, 1, z_score))

as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(sampleTable$condition))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_sig)
for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
   dev.off()
}
par(oma=c(2,2,2,2))
pheatmap(vst_sig_z_scored,
         annotation_col = col_annot,
         cellwidth = 20,
         show_rownames = T,
         scale = "none",
         distfun = "pearson",        
         hclustfun = "ward",
         filename=snakemake@output[["total_heatmap"]], width = 10, height = 20, fontsize =12)

dev.off()
}



col_annot=as.data.frame(as.vector(sampleTable$condition))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst)
z_score <- function(x){(x - mean(x))/sd(x)}
vst_z <- t(apply(vst, 1, z_score))

topVarianceGenes <- head(order(rowVars(assay(vst_dds)), decreasing=T),100)
matrix <- vst_z[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)

for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
   dev.off()
}

par(oma=c(2,2,2,2))


pheatmap(matrix,
         annotation_col = col_annot,
         cellwidth = 20,
         show_rownames = T,
         annotation_names_row=T,
         scale = "none",
         distfun = "pearson",        
         hclustfun = "ward",
         filename=snakemake@output[["total_heatmap"]], width = 10, height = 20, fontsize =12)
dev.off()