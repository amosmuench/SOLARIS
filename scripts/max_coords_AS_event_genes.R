gtf_AS_isoforms=snakemake@params[["gtf_AS_isoforms"]]
library(GenomicFeatures)

txdb <- makeTxDbFromGFF(gtf_AS_isoforms, format="gtf")
genes <- genes(txdb)
genes<-as.data.frame(genes)[,-4]
rownames(genes) <- NULL

write.table(genes, file=snakemake@output[["gtf_genes"]], sep=",")
