library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)
library(topGO)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(ComplexHeatmap)
library(svglite)
subdir<-dirname(snakemake@input[["whippet_mapping_summary"]])
print(subdir)
working_dir=getwd()
print(working_dir)
output_dir=paste(working_dir, "/", subdir, sep="")
print(output_dir)
dir.create(file.path(output_dir, "/GOTerm/"))



########### Mapping Summary ###########
mapped=read.csv(snakemake@input[["whippet_mapping_summary"]], header=T)
mapped[] <- lapply(mapped, gsub, pattern='%', replacement='')
mapped$Mapped_Percent=as.numeric(mapped$Mapped_Percent)
mapped$Multimap_Percent=as.numeric(mapped$Multimap_Percent)
mapped$Novel_Junc_Percent=as.numeric(mapped$Novel_Junc_Percent)
mapped_new=melt(mapped, id.vars = "Run")

mapping_summary=ggplot(mapped_new, aes(factor(Run), y=value, fill = variable)) +
  geom_bar(stat="identity", position = "stack",  aes(factor(Run), y=value, fill = variable)) +
  scale_y_continuous(limits=c(0,100)) +
  xlab("Run") +
  ylab("Total Reads [%]")+
  labs(title = "Read Mapping and Quality",
       subtitle = waiver(),
       caption = waiver(),
       fill = "")+
  scale_fill_manual(labels = c("Mapped Reads", "Multimapped Reads", "Novel Junction Reads"), values=c("forestgreen", "firebrick3", "darkviolet") ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file=snakemake@output[["whippet_mapping_summary_plot"]], plot=mapping_summary, width=10, height=8)

########### AS events frequency ###########
sum_global=read.csv(snakemake@input[["event_table"]], header=T)
################local runs only
#setwd("C:/Users/Amos/Desktop/Master Thesis/Bioinformatic Analysis/scripts")
#sum_global <- read.csv("summary_sig_global.csv",  header=TRUE)
sum_global$comparison<-gsub("_sig_global", "", sum_global$comparison)

########################
#summary(sum_global)

as_freq<- sum_global %>%
  group_by(comparison,Type) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n)*100, digits = 1))

as_event_frequency=ggplot(as_freq, aes(x = comparison, y = freq, label=freq)) +
  geom_bar(aes(color = Type, fill = Type),     stat = "identity" ,position = position_stack() ) +
  geom_text(aes(fill = Type), inherit.aes = TRUE, position = position_stack(vjust=0.5)) +
  ylab(" Frequency of AS Type [%]") +
  xlab("Comparison")+
  ggtitle("Frequency of AS Event Types per Comparison")+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust=0))

ggsave(file=snakemake@output[["whippet_AS_type_frequency_plot"]], plot=as_event_frequency, width=10, height=8)

########### Overlap DAS Genes among conditions ###########

gene_overlap<-sum_global[,c(1,12)]
gene_overlap=gather(gene_overlap, Gene, comparison) %>%
  group_by(Gene, comparison) %>%
  tally %>%
  spread(comparison, n, fill = 0)

gene_overlap=column_to_rownames(gene_overlap, var = "Gene")


gene_overlap[] <- lapply(gene_overlap, function(x) ifelse(x>1, 1, x)) #genes with multiple events are only counted once

gene_overlap_comb=make_comb_mat(gene_overlap,  mode = "intersect")
gene_overlap_comb=gene_overlap_comb[comb_degree(gene_overlap_comb) >= 2]


svg(file=snakemake@output[["whippet_DAS_gene_overlap_among_comparisons_plot"]], height = unit(8, "cm"), width=unit(18, "cm"))
UpSet(gene_overlap_comb, height = unit(8, "cm"), width=unit(18, "cm"),comb_order = order(comb_size(gene_overlap_comb)))
grid.text("Overlap of DAS Genes among Comparisons",x = 0.49, y=0.9, gp=gpar(fontsize=20))
dev.off()

########### Overlap DAS Gene Events among conditions ###########
event_comparison=data.frame(comparison=sum_global$comparison)
event_comparison$event=paste(sum_global$Gene, sum_global$Coord)
event_comparison=gather(event_comparison, event, comparison) %>%
  group_by(event, comparison) %>%
  tally %>%
  spread(comparison, n, fill = 0)
event_comparison=column_to_rownames(event_comparison, var = "event")

event_comparison_comb=make_comb_mat(event_comparison,  mode = "intersect")
event_comparison_comb=event_comparison_comb[comb_degree(event_comparison_comb) >= 2]

svg(file=snakemake@output[["whippet_DAS_event_overlap_among_comparisons_plot"]],  height = unit(8, "cm"), width=unit(18, "cm"))
UpSet(event_comparison_comb, height = unit(8, "cm"), width=unit(18, "cm"),comb_order = order(comb_size(event_comparison_comb)))
grid.text("Overlap of AS events among comparisons",x = 0.49, y=0.9, gp=gpar(fontsize=20))
dev.off()

#################### Prerequisites for GO TERM analysis

#namespace problem tidyr/panther
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

species_gene_goterm_universe_from_panther<-function(taxid)
{#enter taxid, get list of all genes, MISSING RETURN GO TERM TABLE
  detach_package("tidyr")
  library(PANTHER.db)
  panther_species<-availablePthOrganisms(PANTHER.db)
  panther_species<-as.data.frame(panther_species)
  panther_species_name<-subset(panther_species, `UNIPROT Taxon ID`==taxid)# make 3702 passable variable from scheduling script
  pthOrganisms(PANTHER.db)<-panther_species_name$`AnnotationDbi Species`
  go_ids <- keys(PANTHER.db,keytype="GOSLIM_ID")
  cols <- c("ENTREZ","GOSLIM_ID")
  entrez_go_long <-select(PANTHER.db, keys=go_ids, columns=cols,keytype="GOSLIM_ID")
  library(tidyr)
  entrez_go_long=entrez_go_long[!is.na(entrez_go_long$ENTREZ),] #apparently contains na, how?
  library(data.table)
  entrez_go_vector<-setDT(entrez_go_long)[, .(list(unlist(GOSLIM_ID))), ENTREZ]
  entrez_go_vector_list<-as.list(entrez_go_vector$V1)
  names(entrez_go_vector_list)<-entrez_go_vector$ENTREZ
  #str(head(entrez_go_vector_list))
  all_genes<- entrez_go_vector$ENTREZ
  id2GO<-entrez_go_vector_list
  db<-list(id2GO, all_genes)
  names(db)<-c("id2GO","all_genes")
  return(db)

}

#missing: this needs to be taken from snakefile pickle
db<-species_gene_goterm_universe_from_panther("3702")
#################### GO TERM genes AS in single compairisons:

nested_list_aliases2entrez<-function(nested_list_of_genes)
{
  require(mygene)
  genes_event_overlap_entrez=list()
  for(i in names(nested_list_of_genes)){
    df=getGenes(nested_list_of_genes[i], fields="entrezgene")
    df=df[!is.na(df$entrezgene),]
    genes_event_overlap_entrez[[i]]=df$entrezgene
  }
  return(genes_event_overlap_entrez)
}

nested_interesting_genes_and_universe<-function(gene_universe_vector, nested_interesting_genes)
{
  genelist=list()
  for(i in names(nested_interesting_genes)){
    temp=list()
    temp[[i]] <- factor(as.integer (db$all_genes %in% nested_interesting_genes[[i]]))

    names(temp[[i]])<-db$all_genes
    genelist[[i]]<-temp
  }
  return(genelist)
}
go_interesting_genes<-function(nested_interesting_universe, entrez_list, outputname){
  go_results_df=list()
  genes_responsible_for_go_list=list()
  plot_list=list()
  for(i in names(nested_interesting_universe)){
    temp=list()
    temp=nested_interesting_universe[[i]]
    temp_go<-new("topGOdata", ontology="BP", allGenes=temp[[i]], nodeSize=10, annot = annFUN.gene2GO, gene2GO = db$id2GO)
    #graph(test_go)
    test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff=0.05)
    resultFisher <- getSigGroups(temp_go, test.stat)
    resultFisher
    pvalFis <- score(resultFisher)

    #hist(pvalFis, 50, xlab = "p-values")
    allRes <- GenTable(temp_go, classic = resultFisher, topNodes = 100)

    allRes<-as.data.frame(allRes)
    allRes<-subset(allRes, classic<=0.05 & Significant>=5)

    #get the genes responsible for goterm enrichment
    myterms = allRes$GO.ID
    mygenes <- genesInTerm(temp_go, myterms)
    genes_with_go=list()
    for(k in names(mygenes)){
      intermed<-mygenes[[k]]
      intermed<-intermed[intermed %in% entrez_list[[i]]]
      genes_with_go[[k]]<-intermed
    }
    allRes$Significant_Genes<- genes_with_go[match(allRes$GO.ID, names(genes_with_go))]
    allRes$Significant_Genes <- vapply(allRes$Significant_Genes, paste, collapse = ", ", character(1L)) #necessary to write as csv
    allRes$Fold_Enrichment<-allRes$Significant/allRes$Expected
    allRes <-allRes[order(allRes$Fold_Enrichment),]
    allRes$Term <- paste(allRes$GO.ID, allRes$Term, sep=", ")
    allRes$Term <- factor(allRes$Term, levels=rev(allRes$Term))


    if (dim(allRes)[1] != 0) {
      output_second<-i
      if (nchar(output_second)>=70){output_second<-paste(substr(output_second, 1, 70), "\n", substr(output_second, 71,nchar(output_second)))}
      #this should be dynamic
      require(ggplot2)
      go_plot<-ggplot(allRes, aes(x=Term, y=Fold_Enrichment)) +
        stat_summary(geom = "bar", fun = mean, position = "dodge") +
        xlab("Biological process") +
        ylab("Fold Enrichment") +
        ggtitle(paste0("Enriched GO Terms","\n", outputname,"\n",output_second)) +
        scale_y_continuous(breaks = round(seq(0, max(allRes$Fold_Enrichment), by = 0.5), 1)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(
          legend.position='none',
          legend.background=element_rect(),
          plot.title=element_text(angle=0, size=10, face="bold", vjust=1),
          axis.text.x=element_text(angle=0, size=10, hjust=1.10),
          axis.text.y=element_text(angle=0, size=10, vjust=0.5),
          axis.title=element_text(size=12),
          legend.key=element_blank(),
          legend.key.size=unit(1, "cm"),
          legend.text=element_text(size=18),
          title=element_text(size=12)) +
        guides(colour=guide_legend(override.aes=list(size=2.5))) +
        coord_flip()
      plot_list[[i]]<-go_plot
      go_results_df[[i]]<-allRes #should make it into a list of dfs
      #genes_responsible_for_go_list[[i]]<-genes_responsible_for_go
      go_results<-list(plot_list, go_results_df)
      names(go_results)<-c("plot_list", "go_results_df")

    }



  }
  sapply(names(go_results_df),
         function (x) write.table(go_results_df[[x]], file=paste(output_dir, "/GOTerm/",outputname,"_", x, ".txt", sep="") )   )

  for(i in names(plot_list)){
    file_name = paste(output_dir, "/GOTerm/",outputname, "_",i, ".svg", sep="") #_AS_genes_GOTerm
    svglite(file_name, width=10, height=8)
    print(plot_list[[i]])
    dev.off()
  }
  return(go_results)}


as_genes=sum_global[,c(1,12)]


as_genes<-as_genes%>%pivot_wider(id_cols=Gene, values_from="Gene", names_from=comparison)

as_genes_list=list()
for( i in colnames(as_genes)){
  temp=list()
  temp=unnest(as_genes[i], cols=i)
  as_genes_list[[i]]=as.vector(temp[[i]])
}


as_genes_list_entrez<-nested_list_aliases2entrez(as_genes_list)

as_genes_and_universe<-nested_interesting_genes_and_universe(db$all_genes,as_genes_list_entrez)


as_genes_go<-go_interesting_genes(as_genes_and_universe, as_genes_list_entrez, "AS_Genes")







###################################### exact event overlapping among 2 comparisons

event_comparison_comb_2=event_comparison_comb[comb_degree(event_comparison_comb) == 2]
event_comparison_comb_2_names=comb_name(event_comparison_comb_2)

event_list=list()
for (i in event_comparison_comb_2_names){
  event_list[[i]]=extract_comb(event_comparison_comb_2, i)
}

names_matrix=as.data.frame(event_comparison_comb_2)

names_list_1=list()
names_list_2=list()
for( i in colnames(names_matrix)){
  names_list_1[i]=rownames(names_matrix)[(which(names_matrix[i]!=0)[1])]
  names_list_2[i]=rownames(names_matrix)[(which(names_matrix[i]!=0)[2])]
}
#unlist(names_list_1)

#unlist(names_list_2)

comparison_names_overlap=mapply(paste, names_list_1, names_list_2, sep="_AND_", SIMPLIFY=FALSE)

names(event_list)<-comparison_names_overlap
temp=list()
genes_event_overlap=list()
for(i in names(event_list)){
  temp[i]= event_list[i]
  for(j in temp[i]){
    splitted=list()
    splitted=lapply(strsplit(j," "), function(x) x[1])
    splitted=unique(splitted)
    genes_event_overlap[[i]]=unlist(splitted)
  }
}


genes_event_overlap_entrez<-nested_list_aliases2entrez(genes_event_overlap)

event_overlap_universe<-nested_interesting_genes_and_universe(db$all_genes, genes_event_overlap_entrez)


event_overlap_go<-go_interesting_genes(event_overlap_universe, genes_event_overlap_entrez, "AS_event_overlap")



################################### gene overlap among comparisons

gene_overlap_comb=make_comb_mat(gene_overlap,  mode = "intersect")
gene_overlap_comb=gene_overlap_comb[comb_degree(gene_overlap_comb) >= 2]

gene_overlap_comb_2=gene_overlap_comb[comb_degree(gene_overlap_comb) == 2]

gene_overlap_comb_2_names=comb_name(gene_overlap_comb_2)
gene_overlap_list=list()
for (i in gene_overlap_comb_2_names){
  gene_overlap_list[[i]]=extract_comb(gene_overlap_comb_2, i)
}

names_matrix=as.data.frame(gene_overlap_comb_2)

names_list_1=list()
names_list_2=list()
for( i in colnames(names_matrix)){
  names_list_1[i]=rownames(names_matrix)[(which(names_matrix[i]!=0)[1])]
  names_list_2[i]=rownames(names_matrix)[(which(names_matrix[i]!=0)[2])]
}
#unlist(names_list_1)

#unlist(names_list_2)

comparison_names_overlap=mapply(paste, names_list_1, names_list_2, sep="_AND_", SIMPLIFY=FALSE)


names(gene_overlap_list)<-comparison_names_overlap




gene_overlap_entrez<-nested_list_aliases2entrez(gene_overlap_list) #this is somehow broken

genes_overlap_universe<-nested_interesting_genes_and_universe(db$all_genes, gene_overlap_entrez)


genes_overlap_go<-go_interesting_genes(genes_overlap_universe, gene_overlap_entrez, "AS_gene_overlap")
