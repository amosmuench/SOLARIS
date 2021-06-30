library(DESeq2)
library(rhdf5)
library("tximport")
library(GenomicFeatures)
library(mygene)
library(ggplot2)
library(svglite)
library(pheatmap)


abundance_files<-snakemake@input[["abundance_h5"]]
dge_names<-snakemake@params[["dge_names"]]
print(abundance_files)
print(dge_names)

dge_names=gsub(" ", "", dge_names, fixed = TRUE)
print(dge_names)
dge_names=make.names(dge_names)
gtf<-snakemake@params[["gff_dge"]]
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
#t<-keys(txdb, keytype = "TXNAME")
#g<-keys(txdb, keytype = "GENEID")
#t<-as.data.frame(t)
#head(t)
#g<-as.data.frame(g)
#t$gene<-sub("\\.\\d+$", "", t$t)

#colnames(t)<-c("TXNAME", "GENEID")
print("revised tx2gene")
head(tx2gene)

imported_ab<-tximport(files, type="kallisto", tx2gene = tx2gene)

sampleTable<-data.frame(condition=factor(dge_names))
print(sampleTable)

#rownames(sampleTable)<-colnames(imported_ab$counts)

rownames(sampleTable)=make.names(colnames(imported_ab$counts), unique=TRUE)






dds<-DESeqDataSetFromTximport(imported_ab, sampleTable, ~condition)
vst_dds <- vst(dds, blind=TRUE, fitType = "local")
vst=assay(vst_dds)



#make names and dirs for plots
cwd=getwd()
study_name=snakemake@params[["study_name"]]
dir.create(file.path(cwd, study_name, "DGE_DTE/genes_of_interest"))
outpath=file.path(cwd, study_name, "DGE_DTE/genes_of_interest")
comp_name=paste(dge_names, collapse = '')

goi=snakemake@params[["goi"]]

print(goi)
goi=gsub("\\[|\\]|\\'", "",goi) #|\\'|\\,",

goi=sapply(goi, strsplit, ",")
goi=goi[[1]]
print(goi)
print(dds)
print(goi[goi %in% rownames(dds)])
goi=goi[goi %in% rownames(dds)]
goi_symbols=getGenes(goi, field="symbol")
print("mistake?")
print(goi)
print(goi_symbols)
sample_goi_plotcounts=list()
#retrieve normalization factor manually
#co=counts(dds)
#co<-as.data.frame(co)
#head(co)
#cou2=co
#head(cou2)
#print("after")
#print(colnames(cou2))
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#gm_gene=vector()
#for(i in 1:nrow(cou2)) {
#    row <- cou2[i,]
#    gm_gene<-c(gm_gene,gm_mean(row))   
#}
#cou2$pseudo_ref=gm_gene
#head(cou2)
#norm_factor=vector()
#for(i in dge_names)
#{
#cou2[[paste0("ratio_",i,"_2ref")]]=cou2[[toString(i)]]/cou2[["pseudo_ref"]]
#norm_factor=c(norm_factor, median(cou2[[paste0("ratio_",i,"_2ref")]]))
#}

#print(norm_factor)
for(i in goi)
{
print(i)
 sample_goi_plotcounts[[i]]=plotCounts(dds, gene=i, returnData = T, transform=F, normalized=T)
}


goi_plots=list()
goi_symbols=goi_symbols[,c("query", "symbol")]
print(goi_symbols)
goi_symbols=unique(goi_symbols)
print("here")
print(list(names(sample_goi_plotcounts)))
goi_symbols=goi_symbols[match(names(sample_goi_plotcounts), goi_symbols$query),]

transcript_symbol=c(list(names(sample_goi_plotcounts)), list(unique(goi_symbols$symbol)))
print(transcript_symbol)
for(i in seq_along(transcript_symbol[[1]]))
{
  goi_plots[[paste0(transcript_symbol[[2]][[i]],"_", transcript_symbol[[1]][[i]])]]=
    ggplot(data=sample_goi_plotcounts[[i]], aes(x=condition, y=count))+
    geom_boxplot()+
    ggtitle(paste0(transcript_symbol[[2]][[i]], "\n", transcript_symbol[[1]][[i]] ))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #print(goi_plots[[paste0(transcript_symbol[[2]][[i]],"_", transcript_symbol[[1]][[i]])]])
}

for(i in names(goi_plots))
{
  ggsave(file=file.path(outpath, paste0("_", i, "_count_plot.png")), plot=goi_plots[[i]])
}
write.table(goi, snakemake@output[[1]])



#heatmaps goi

dds<-DESeqDataSetFromTximport(imported_ab, sampleTable, ~condition)
vst_dds <- vst(dds, blind=TRUE, fitType = "local")
vst=assay(vst_dds)

contrast_combinations=combn(levels(sampleTable$condition), 2, FUN=list)
comb_list_goi=list()
for(i in contrast_combinations)
{
print(i[1])
print(i[2])
print(head(vst))

sub_mat <- matrix(c(i[1], i[2]), ncol = 1)

vst_goi=vst[unlist(goi), ]
comb_list_goi[[paste0(i[1],"_vs_", i[2])]]<- vst_goi[, sub_mat[, 1]]


comb_name=paste0(i[1],"_vs_", i[2])

print(comb_list_goi[[comb_name]])

#genenames_split=rownames(comb_list_goi[[comb_name]])
 
query_list=rownames(comb_list_goi[[comb_name]])
print(query_list)
gene_symbols=getGenes(query_list, field="symbol")
print("here")
print(gene_symbols)
gene_symbols<- ifelse(is.na(gene_symbols$symbol), gene_symbols$query, gene_symbols$symbol)

goi_symbols=goi_symbols[,c("query", "symbol")]
print(goi_symbols)
goi_symbols=unique(goi_symbols)
print("here")
goi_symbols=goi_symbols[match(goi, goi_symbols$query),]
print(goi_symbols)
print(length(goi_symbols$symbol))
print(length(unique(goi_symbols$symbol)))

print(rownames(comb_list_goi[[comb_name]]))
print(length(rownames(comb_list_goi[[comb_name]])))
print(goi_symbols$symbol)


rownames(comb_list_goi[[comb_name]])<-goi_symbols$symbol

z_score <- function(x){(x - mean(x))/sd(x)}
as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(c(i[1], i[2]) ))
colnames(col_annot)="condition"
print("cond")
rownames(col_annot)<-colnames(comb_list_goi[[comb_name]])
vst_goi_z_scored <- t(apply(comb_list_goi[[comb_name]], 1, z_score))

#t(scale(t(log2matrix)))
#as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(c(i[1], i[2])))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_goi_z_scored)


par(oma=c(2,2,2,2))
print(file.path(outpath, paste0("heatmap_goi_", i[1], "_vs_", i[2], ".png")))
pheatmap(vst_goi_z_scored,
         annotation_col = col_annot,
         cellwidth = 20,
         show_rownames = T,
         annotation_names_row=T,
         scale = "none",
         distfun = "pearson",        
         hclustfun = "ward",
         filename=file.path(outpath, paste0("heatmap_goi_", i[1], "_vs_", i[2], ".pdf")), width = 10, height = 20, fontsize =12)
dev.off()



l2fc_goi=log2(comb_list_goi[[comb_name]][,c(i[1])])-log2(comb_list_goi[[comb_name]][,c(i[2])])
#l2fc_=log2(comb_list_goi[[comb_name]][,c(i[2])]/comb_list_goi[[comb_name]][,c(i[2])]) #wrong
print(l2fc_goi)
print(rownames(comb_list_goi[[comb_name]]))
l2fc_table=as.data.frame(as.list(l2fc_goi) )
write.table(l2fc_table, file=file.path(outpath, paste0(comb_name, "goi_l2fc.tsv")), sep="\t")

######global
comb_list=list()
comb_list[[paste0(i[1],"_vs_", i[2])]]<- vst[, sub_mat[, 1]]



z_score <- function(x){(x - mean(x))/sd(x)}
as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(c(i[1], i[2]) ))
colnames(col_annot)="condition"
print("cond")
rownames(col_annot)<-colnames(comb_list[[comb_name]])
vst_z_scored <- t(apply(comb_list[[comb_name]], 1, z_score))
print(dim(vst_z_scored))
vst_z_scored=na.omit(vst_z_scored)
print(dim(vst_z_scored))
#t(scale(t(log2matrix)))
#as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(c(i[1], i[2])))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_z_scored)



par(oma=c(2,2,2,2))
print(file.path(outpath, paste0("heatmap_all_", i[1], "_vs_", i[2], ".pdf")))
#pheatmap(vst_z_scored,
#         annotation_col = col_annot,
#         cellwidth = 20,
#         show_rownames = T,
#         annotation_names_row=T,
#         scale = "none",
#         distfun = "pearson",        
#         hclustfun = "ward",
#         filename=file.path(outpath, paste0("heatmap_all_", i[1], "_vs_", i[2], ".pdf")), width = 10, height = 20, fontsize =12)
#dev.off()


l2fc_all=log2(comb_list[[comb_name]][,c(i[1])])-log2(comb_list[[comb_name]][,c(i[2])])
#l2fc_=log2(comb_list_goi[[comb_name]][,c(i[2])]/comb_list_goi[[comb_name]][,c(i[2])]) #wrong

#l2fc_table_all=as.data.frame(as.list(l2fc_all) )
write.table(l2fc_all, file=file.path(outpath, paste0(comb_name, "_all_l2fc.tsv")), sep="\t")

}

vst_goi=vst[unlist(goi), ]
z_score <- function(x){(x - mean(x))/sd(x)}
as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(sampleTable$condition))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_goi)
vst_goi_z_scored <- t(apply(vst_goi, 1, z_score))

#t(scale(t(log2matrix)))
as.vector(sampleTable$condition)
col_annot=as.data.frame(as.vector(sampleTable$condition))
colnames(col_annot)="condition"
rownames(col_annot)<-colnames(vst_goi_z_scored)


par(oma=c(2,2,2,2))


pheatmap(vst_goi_z_scored,
         annotation_col = col_annot,
         cellwidth = 20,
         show_rownames = T,
         annotation_names_row=T,
         scale = "none",
         distfun = "pearson",        
         hclustfun = "ward",
         filename=snakemake@output[["goi_heatmap_all_cond"]], width = 10, height = 20, fontsize =12)
dev.off()

