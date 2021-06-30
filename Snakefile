import pickle
import os
from collections import OrderedDict

cwd =os.getcwd()
def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


				
#snakemake functions only take 1 wildcard argument hence duplication
def replicates_whippet_A(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["A"]
     
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin)  # + ","not tested, perhaps needs extra ,
       sample_begin=",".join(sample_begin)+","
       return sample_begin
     else:
       sample_list=list()
       for i in range(len(sample_begin)):
         
         each=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       sample_list=",".join(sample_list)+","
       print("sample_list_A-re", sample_list)
       return sample_list        

def replicates_whippet_B(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["B"]
     
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin) 
       sample_begin=",".join(sample_begin)+","
       return sample_begin
     else:
       sample_list=list()
       for i in range(len(sample_begin)):
         
         each=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       sample_list=",".join(sample_list)+","
       
       return sample_list 

 
def input_whippet_delta_A(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["A"]
     
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin)
       return sample_begin
     else:
       sample_list=list()
       
       for i in range(len(sample_begin)):
         each=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       
       return sample_list        
 
def input_whippet_delta_B(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["B"]
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin)
       return sample_begin
     else:
       sample_list=list()
       print(sample_begin)
       for i in range(len(sample_begin)):
         each=expand("{study}/AS/{sample}.psi.gz", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       return sample_list 

    
def input_bam_A(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["A"]
     
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/BAM_sorted/{sample}.bam", study=studies, sample=sample_begin)
       return sample_begin
     else:
       sample_list=list()
       
       for i in range(len(sample_begin)):
         each=expand("{study}/AS/BAM_sorted/{sample}.bam", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       
       return sample_list 

def input_bam_B(wildcards):
     sample_begin= comparison_dict[wildcards.comparison_whippet]["B"]
     
     if isinstance(sample_begin, str):
       sample_begin=expand("{study}/AS/BAM_sorted/{sample}.bam", study=studies, sample=sample_begin)
       return sample_begin
     else:
       sample_list=list()
       
       for i in range(len(sample_begin)):
         each=expand("{study}/AS/BAM_sorted/{sample}.bam", study=studies, sample=sample_begin[i])
         each=''.join(each)
         sample_list.append(each)
       print("bam_list_B", sample_list)
       return sample_list 

items = loadall("study_total.pkl")
directory= cwd

studies=next(items)
tax_id=next(items)
read_files=next(items)
paired=next(items)
replicates=next(items)
complete_list=next(items)
comparison_dict=next(items)
comparison_dict_v2=next(items)
comp_A_names=next(items)
sample_acc_list_A=next(items)
run_accession_list_A=next(items)
run_list_A=next(items)
comp_B_names=next(items)
sample_acc_list_B=next(items)
run_accession_list_B=next(items)
run_list_B=next(items)

contrast_names_whippet = [i + "_vs_"+ j for i, j in zip(comp_A_names, comp_B_names)]



try:
  dge_runs=list(OrderedDict.fromkeys(run_accession_list_A))
except TypeError:
  dge_runs = [tuple(l) for l in run_accession_list_A]
dge_names= list(OrderedDict.fromkeys(comp_A_names)) #this is false as soon as number of comparisons is uneven, better: runacc to list and as pickle, order is saved



###############input function for whippet
groupA=dict(zip(comp_A_names, sample_acc_list_A))

groupB=dict(zip(comp_B_names, sample_acc_list_B))

comparison_whippet_wc=comparison_dict.keys()


configfile: "/home/amosmuench/Pipeline17092020/config.yaml" #make this dependent on cwd
fasta= config["fasta"]
gtf_AS_isoforms= config["gtf_AS_isoforms"]
julia_dir=config["julia_dir"]
whippet_dir=config["whippet_dir"]
gff_dge=config["gff_dge"]


whippet_out_psi=expand("{study}/AS/{read}.psi.gz", study=studies, read=read_files)
print("whippet_out_psi", whippet_out_psi)

if any(isinstance(el, list) for el in read_files):

  rep=len(read_files[0])

  read_files = [item for sublist in read_files for item in sublist]

  dge_names=[i for j in dge_names for i in (j,)*rep]



  

if paired:
    if config["gene_ids_of_interest"]:
        rule all:
            input:
                whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
                whippet_out_psi=expand("{study}/AS/{read}.psi.gz", study=studies, read=read_files),
                whippet_out_gene_tpm=expand("{study}/AS/{read}.gene.tpm.gz", study=studies, read=read_files),
                whippet_out_isoform_tpm=expand("{study}/AS/{read}.isoform.tpm.gz", study=studies, read=read_files),
                whippet_out_jnc=expand("{study}/AS/{read}.jnc.gz", study=studies, read=read_files),
                whippet_out_map=expand("{study}/AS/{read}.map.gz", study=studies, read=read_files),
                whippet_out_sam=expand("{study}/AS/{read}.sam", study=studies, read=read_files),
                multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies),
                comparison_whippet_wc=expand("{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz", study=studies, comparison_whippet=comparison_whippet_wc),
                kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
                whippet_out_sorted=expand("{study}/AS/BAM_sorted/{read}.bam"   , study=studies, read = read_files),
                whippet_out_sorted_indexed=expand("{study}/AS/BAM_sorted/{read}.bam.bai", study=studies, read = read_files),
                abundance_h5=expand("{study}/DGE_DTE/quant/{read}/abundance.h5", study=studies, read=read_files),
                abundance_tsv=expand("{study}/DGE_DTE/quant/{read}/abundance.tsv", study=studies, read=read_files),
                run_info=expand("{study}/DGE_DTE/quant/{read}/run_info.json", study=studies, read=read_files),
                #star_out=expand("{study}/Alignment/{read}/Aligned.out.sam",study=studies, read=read_files),
                #star_index=expand("{study}/Alignment/star_index/SAindex", study=studies),
                diff_unzipped=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                diff_sig_global=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff",study=studies,  comparison_whippet=comparison_whippet_wc),
                diff_sig_goi=expand("{study}/AS/DeltaPSI/genes_of_interest/{comparison_whippet}_sig_goi.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                whippet_mapping_dc=expand("{study}/AS/{read}.map", study=studies, read=read_files),
                whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies),
                event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies),
                whippet_mapping_summary_plot=expand("{study}/AS/Summary/Whippet_mapping_summary_plot.svg", study=studies),
                whippet_AS_type_frequency_plot=expand("{study}/AS/Summary/Whippet_Frequency_AS_Event_Types_plot.svg", study=studies),
                whippet_DAS_gene_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_gene_overlap_among_comparisons_plot.svg", study=studies),
                whippet_DAS_event_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_event_overlap_among_comparisons_plot.svg", study=studies),
                bam_tsv=expand("{study}/AS/DeltaPSI/{comparison_whippet}_bam.tsv",study=studies, comparison_whippet=comparison_whippet_wc),
                mock_plot=expand("{study}/AS/sashimi_plots/{comparison_whippet}.txt",study=studies, comparison_whippet=comparison_whippet_wc),
                mock_plot_detail=expand("{study}/AS/sashimi_plots/{comparison_whippet}_detail.txt",study=studies, comparison_whippet=comparison_whippet_wc),
                total_heatmap=expand("{study}/DGE_DTE/Summary/total_heatmap.png", study=studies),
                pca=expand("{study}/DGE_DTE/Summary/pca.png", study=studies),
                mock_gene_plot=expand("{study}/DGE_DTE/geneplots.txt", study=studies),
                goi_heatmap_all_cond=expand("{study}/DGE_DTE/Summary/goi_all_cond_heatmap.png", study=studies)
    else:
        rule all:
            input:
                whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
                whippet_out_psi=expand("{study}/AS/{read}.psi.gz", study=studies, read=read_files),
                whippet_out_gene_tpm=expand("{study}/AS/{read}.gene.tpm.gz", study=studies, read=read_files),
                whippet_out_isoform_tpm=expand("{study}/AS/{read}.isoform.tpm.gz", study=studies, read=read_files),
                whippet_out_jnc=expand("{study}/AS/{read}.jnc.gz", study=studies, read=read_files),
                whippet_out_map=expand("{study}/AS/{read}.map.gz", study=studies, read=read_files),
                whippet_out_sam=expand("{study}/AS/{read}.sam", study=studies, read=read_files),
                multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies),
                comparison_whippet_wc=expand("{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz", study=studies, comparison_whippet=comparison_whippet_wc),
                kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
                whippet_out_sorted=expand("{study}/AS/BAM_sorted/{read}.bam"   , study=studies, read = read_files),
                whippet_out_sorted_indexed=expand("{study}/AS/BAM_sorted/{read}.bam.bai", study=studies, read = read_files),
                abundance_h5=expand("{study}/DGE_DTE/quant/{read}/abundance.h5", study=studies, read=read_files),
                abundance_tsv=expand("{study}/DGE_DTE/quant/{read}/abundance.tsv", study=studies, read=read_files),
                run_info=expand("{study}/DGE_DTE/quant/{read}/run_info.json", study=studies, read=read_files),
                #star_out=expand("{study}/Alignment/{read}/Aligned.out.sam",study=studies, read=read_files),
                #star_index=expand("{study}/Alignment/star_index/SAindex", study=studies),
                diff_unzipped=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                diff_sig_global=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                whippet_mapping_dc=expand("{study}/AS/{read}.map", study=studies, read=read_files),
                whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies),
                event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies),
                whippet_mapping_summary_plot=expand("{study}/AS/Summary/Whippet_mapping_summary_plot.svg", study=studies),
                whippet_AS_type_frequency_plot=expand("{study}/AS/Summary/Whippet_Frequency_AS_Event_Types_plot.svg", study=studies),
                whippet_DAS_gene_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_gene_overlap_among_comparisons_plot.svg", study=studies),
                whippet_DAS_event_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_event_overlap_among_comparisons_plot.svg", study=studies),
                total_heatmap=expand("{study}/DGE_DTE/Summary/total_heatmap.png", study=studies),
                pca=expand("{study}/DGE_DTE/Summary/pca.png", study=studies)










else:
    if config["gene_ids_of_interest"]:
        rule all:
            input:
                whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
                whippet_out_psi=expand("{study}/AS/{read}.psi.gz", study=studies, read=read_files),
                whippet_out_gene_tpm=expand("{study}/AS/{read}.gene.tpm.gz", study=studies, read=read_files),
                whippet_out_isoform_tpm=expand("{study}/AS/{read}.isoform.tpm.gz", study=studies, read=read_files),
                whippet_out_jnc=expand("{study}/AS/{read}.jnc.gz", study=studies, read=read_files),
                whippet_out_map=expand("{study}/AS/{read}.map.gz", study=studies, read=read_files),
                whippet_out_sam=expand("{study}/AS/{read}.sam", study=studies, read=read_files),
                multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies),
                comparison_whippet_wc=expand("{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz", study=studies, comparison_whippet=comparison_whippet_wc),
                kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
                whippet_out_sorted=expand("{study}/AS/BAM_sorted/{read}.bam"   , study=studies, read = read_files),
                whippet_out_sorted_indexed=expand("{study}/AS/BAM_sorted/{read}.bam.bai", study=studies, read = read_files),
                abundance_h5=expand("{study}/DGE_DTE/quant/{read}/abundance.h5", study=studies, read=read_files),
                abundance_tsv=expand("{study}/DGE_DTE/quant/{read}/abundance.tsv", study=studies, read=read_files),
                run_info=expand("{study}/DGE_DTE/quant/{read}/run_info.json", study=studies, read=read_files),
                #star_out=expand("{study}/Alignment/{read}/Aligned.out.sam",study=studies, read=read_files),
                #star_index=expand("{study}/Alignment/star_index/SAindex", study=studies),
                diff_unzipped=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}.diff",study=studies, read=read_files, comparison_whippet=comparison_whippet_wc),
                diff_sig_global=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                diff_sig_goi=expand("{study}/AS/DeltaPSI/genes_of_interest/{comparison_whippet}_sig_goi.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                whippet_mapping_dc=expand("{study}/AS/{read}.map", study=studies, read=read_files),
                whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies),
                event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies),
                whippet_mapping_summary_plot=expand("{study}/AS/Summary/Whippet_mapping_summary_plot.svg", study=studies),
                whippet_AS_type_frequency_plot=expand("{study}/AS/Summary/Whippet_Frequency_AS_Event_Types_plot.svg", study=studies),
                whippet_DAS_gene_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_gene_overlap_among_comparisons_plot.svg", study=studies),
                whippet_DAS_event_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_event_overlap_among_comparisons_plot.svg", study=studies),
                bam_tsv=expand("{study}/AS/DeltaPSI/{comparison_whippet}_bam.tsv",study=studies, comparison_whippet=comparison_whippet_wc),
                mock_plot=expand("{study}/AS/sashimi_plots/{comparison_whippet}.txt",study=studies, comparison_whippet=comparison_whippet_wc),
                gtf_genes=expand("{study}/AS/DeltaPSI/genes_of_interest/gtf_genes_coordinates.txt", study=studies),
                mock_plot_detail=expand("{study}/AS/sashimi_plots/{comparison_whippet}_detail.txt",study=studies, comparison_whippet=comparison_whippet_wc),
                total_heatmap=expand("{study}/DGE_DTE/Summary/total_heatmap.png", study=studies),
                pca=expand("{study}/DGE_DTE/Summary/pca.png", study=studies),
                mock_gene_plot=expand("{study}/DGE_DTE/geneplots.txt", study=studies),
                goi_heatmap_all_cond=expand("{study}/DGE_DTE/Summary/goi_all_cond_heatmap.png", study=studies)
    else:
        rule all:
            input:
                whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
                whippet_out_psi=expand("{study}/AS/{read}.psi.gz", study=studies, read=read_files),
                whippet_out_gene_tpm=expand("{study}/AS/{read}.gene.tpm.gz", study=studies, read=read_files),
                whippet_out_isoform_tpm=expand("{study}/AS/{read}.isoform.tpm.gz", study=studies, read=read_files),
                whippet_out_jnc=expand("{study}/AS/{read}.jnc.gz", study=studies, read=read_files),
                whippet_out_map=expand("{study}/AS/{read}.map.gz", study=studies, read=read_files),
                whippet_out_sam=expand("{study}/AS/{read}.sam", study=studies, read=read_files),
                multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies),
                comparison_whippet_wc=expand("{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz", study=studies, comparison_whippet=comparison_whippet_wc),
                kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
                whippet_out_sorted=expand("{study}/AS/BAM_sorted/{read}.bam"   , study=studies, read = read_files),
                whippet_out_sorted_indexed=expand("{study}/AS/BAM_sorted/{read}.bam.bai", study=studies, read = read_files),
                abundance_h5=expand("{study}/DGE_DTE/quant/{read}/abundance.h5", study=studies, read=read_files),
                abundance_tsv=expand("{study}/DGE_DTE/quant/{read}/abundance.tsv", study=studies, read=read_files),
                run_info=expand("{study}/DGE_DTE/quant/{read}/run_info.json", study=studies, read=read_files),
                diff_unzipped=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                diff_sig_global=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff",study=studies, comparison_whippet=comparison_whippet_wc),
                whippet_mapping_dc=expand("{study}/AS/{read}.map", study=studies, read=read_files),
                whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies),
                event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies),
                whippet_mapping_summary_plot=expand("{study}/AS/Summary/Whippet_mapping_summary_plot.svg", study=studies),
                whippet_AS_type_frequency_plot=expand("{study}/AS/Summary/Whippet_Frequency_AS_Event_Types_plot.svg", study=studies),
                whippet_DAS_gene_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_gene_overlap_among_comparisons_plot.svg", study=studies),
                whippet_DAS_event_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_event_overlap_among_comparisons_plot.svg", study=studies),
                total_heatmap=expand("{study}/DGE_DTE/Summary/total_heatmap.png", study=studies)
if paired:
    rule ena_download_paired:
        output:
            expand("{study}/reads/{study}/{read}/{read}_1.fastq.gz", study=studies, read=read_files),
            expand("{study}/reads/{study}/{read}/{read}_2.fastq.gz", study=studies, read=read_files)
        conda:
            "./envs/ena_download.yaml"
        shell:
            f"python ./scripts/enaBrowserTools-1.6/python3/enaGroupGet.py {studies} -f fastq -d {studies}/reads/ echo $"
else:
    rule ena_download_single:
        output:
            expand("{study}/reads/{study}/{read}/{read}.fastq.gz", study=studies, read=read_files)
        conda:
            "./envs/ena_download.yaml"
        shell:
            f"python ./scripts/enaBrowserTools-1.6/python3/enaGroupGet.py {studies} -f fastq -d {studies}/reads/"

if paired:
    rule fastqc_paired:
        input:
            fwd=expand("{study}/reads/{study}/{read}/{read}_1.fastq.gz", study=studies, read=read_files),
            rev=expand("{study}/reads/{study}/{read}/{read}_2.fastq.gz", study=studies, read=read_files)
        output:
            fastqc_fwd_report=expand("{study}/QC/{read}_1_fastqc.html", study=studies, read=read_files),
            fastqc_fwd_zip=expand("{study}/QC/{read}_1_fastqc.zip", study=studies, read=read_files),
            fastqc_rev_report=expand("{study}/QC/{read}_2_fastqc.html", study=studies, read=read_files),
            fastqc_rev_zip=expand("{study}/QC/{read}_2_fastqc.zip", study=studies, read=read_files)
        conda:
            "./envs/fastqc.yaml"
        params:
            fastqc_outdir=expand("{study}/QC/", study=studies)
        threads: 4
        shell:
            "fastqc -t 4 {input.fwd} {input.rev} --outdir={params.fastqc_outdir}"

else:
    rule fastqc_single:
        input:
            expand("{study}/reads/{study}/{read}/{read}.fastq.gz", study=studies, read=read_files),
        output:
            fastqc_report=expand("{study}/QC/{read}_fastqc.html", study=studies, read=read_files),
            fastqc_zip=expand("{study}/QC/{read}_fastqc.zip", study=studies, read=read_files)
        conda:
            "./envs/fastqc.yaml"
        params:
            fastqc_outdir=expand("{study}/QC/", study=studies)
        threads: 4
        shell:
            "fastqc -t 4 {input} --outdir={params.fastqc_outdir}"



if paired:
    rule multiqc_paired:
        input:
            fastqc_fwd_report=expand("{study}/QC/{read}_1_fastqc.html", study=studies, read=read_files),
            fastqc_fwd_zip=expand("{study}/QC/{read}_1_fastqc.zip", study=studies, read=read_files),
            fastqc_rev_report=expand("{study}/QC/{read}_2_fastqc.html", study=studies, read=read_files),
            fastqc_rev_zip=expand("{study}/QC/{read}_2_fastqc.zip", study=studies, read=read_files)
        output:
            multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies)
        conda:
            "./envs/multiqc.yaml"
        params:
            multiqc_outdir=expand("{study}/QC_aggregated/", study=studies),
            multiqc_inputdir=expand("{study}/QC/", study=studies)
        shell:
            "multiqc {params.multiqc_inputdir} -o {params.multiqc_outdir}"

else:
    rule multiqc_single:
        input:
            fastqc_report=expand("{study}/QC/{read}_fastqc.html", study=studies, read=read_files),
            fastqc_zip=expand("{study}/QC/{read}_fastqc.zip", study=studies, read=read_files)
        output:
            multiqc_report=expand("{study}/QC_aggregated/multiqc_report.html", study=studies)
        conda:
            "./envs/multiqc.yaml"
        params:
            multiqc_outdir=expand("{study}/QC_aggregated/", study=studies),
            multiqc_inputdir=expand("{study}/QC/", study=studies)
        shell:
            "multiqc {params.multiqc_inputdir} -o {params.multiqc_outdir}"

rule whippet_index:
    priority: 1
    output:
        whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies) #add to rule all to trigger as it is not connected to any of the other rules and only runs once
    params:
        fasta= config["fasta"],
        gtf_AS_isoforms= config["gtf_AS_isoforms"],
        julia_dir=config["julia_dir"],
        whippet_dir=config["whippet_dir"],
        whippet_index_outdir=expand("{study}/AS/index/AS_index.jls", study=studies)
    shell:
        "{params.julia_dir} {params.whippet_dir}/whippet-index.jl --fasta {params.fasta} --gtf {params.gtf_AS_isoforms} --index {params.whippet_index_outdir}"

if paired:
    rule whippet_quant_paired:
        input:
            whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
            fwd="{study}/reads/{study}/{read}/{read}_1.fastq.gz",
            rev="{study}/reads/{study}/{read}/{read}_2.fastq.gz"
        output:
            whippet_out_psi="{study}/AS/{read}.psi.gz",
            whippet_out_gene_tpm="{study}/AS/{read}.gene.tpm.gz",
            whippet_out_isoform_tpm="{study}/AS/{read}.isoform.tpm.gz",
            whippet_out_jnc="{study}/AS/{read}.jnc.gz",
            whippet_out_map="{study}/AS/{read}.map.gz",
            whippet_out_sam="{study}/AS/{read}.sam"
        params:
            whippet_out_name="{study}/AS/{read}",
            whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
            julia_dir=config["julia_dir"],
            whippet_dir=config["whippet_dir"],
            whippet_index_outdir=expand("{study}/AS/index/AS_index.jls", study=studies),
        shell:
            "{params.julia_dir} {params.whippet_dir}/whippet-quant.jl {input.fwd} {input.rev} -o {params.whippet_out_name} -x {params.whippet_index} --sam > {output.whippet_out_sam}"

else:
    rule whippet_quant_single:
        input:
            whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
            reads="{study}/reads/{study}/{read}/{read}.fastq.gz"
        output:
            whippet_out_psi="{study}/AS/{read}.psi.gz",
            whippet_out_gene_tpm="{study}/AS/{read}.gene.tpm.gz",
            whippet_out_isoform_tpm="{study}/AS/{read}.isoform.tpm.gz",
            whippet_out_jnc="{study}/AS/{read}.jnc.gz",
            whippet_out_map="{study}/AS/{read}.map.gz",
            whippet_out_sam="{study}/AS/{read}.sam"
        params:
            whippet_out_name="{study}/AS/{read}",
            whippet_index=expand("{study}/AS/index/AS_index.jls", study=studies),
            julia_dir=config["julia_dir"],
            whippet_dir=config["whippet_dir"],
            whippet_index_outdir=expand("{study}/AS/index/", study=studies)
        shell:
            "{params.julia_dir} {params.whippet_dir}/whippet-quant.jl {input.reads} -o {params.whippet_out_name} -x {params.whippet_index} --sam > {output.whippet_out_sam}"


rule whippet_delta:
    input:
        input_whippet_delta_A,
        input_whippet_delta_B
        #lambda wildcards : expand("{study}/AS/{sample}.psi.gz", sample= comparison_dict[wildcards.comparison_whippet]["A"], study=studies),
        #lambda wildcards : expand("{study}/AS/{sample}.psi.gz", sample= comparison_dict[wildcards.comparison_whippet]["B"], study=studies)
    output:
        "{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz"
    params:
        condition_A = replicates_whippet_A,
        condition_B = replicates_whippet_B,
        julia_dir=config["julia_dir"],
        whippet_dir=config["whippet_dir"],
        output = lambda wildcards : expand("{study}/AS/DeltaPSI/" + wildcards.comparison_whippet, study=studies)
    shell:
        "{params.julia_dir} {params.whippet_dir}/whippet-delta.jl -a {params.condition_A} -b {params.condition_B} -o {params.output}"

rule bam_sort:
    input:
        whippet_out_sam="{study}/AS/{read}.sam"
    output:
        whippet_out_sorted="{study}/AS/BAM_sorted/{read}.bam"
    conda:
        "./envs/samtools.yaml"
    params:
        threads=8
    shell:
        "samtools sort -@ {params.threads} -o {output.whippet_out_sorted} {input} "

rule bam_index:
    input:
        whippet_out_sorted="{study}/AS/BAM_sorted/{read}.bam"
    output:
        whippet_out_sorted_indexed="{study}/AS/BAM_sorted/{read}.bam.bai"
    conda:
        "./envs/samtools.yaml"
    params:
        threads=8
    shell:
        "samtools index -@ {params.threads} {input[0]}"

rule kallisto_index:
    output:
        kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies)
    params:
        cdna_fasta= config["cdna_fasta"],
        outdir= expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies)
    conda:
        "./envs/kallisto.yaml"
    threads: 8
    priority: 1
    log:
        expand("{study}/logs/kallisto_index.log", study=studies)
    shell:
        "kallisto index --index {params.outdir} {params.cdna_fasta} 2> {log}"

if paired:
    rule kallisto_quant_paired:
        input:
            fwd="{study}/reads/{study}/{read}/{read}_1.fastq.gz",
            rev="{study}/reads/{study}/{read}/{read}_2.fastq.gz"
        output:
            abundance_h5="{study}/DGE_DTE/quant/{read}/abundance.h5",
            abundance_tsv="{study}/DGE_DTE/quant/{read}/abundance.tsv",
            run_info="{study}/DGE_DTE/quant/{read}/run_info.json"
        params:
            kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
            outdir="{study}/DGE_DTE/quant/{read}/",
            threads=8
        conda:
            "./envs/kallisto.yaml"
        log:
            "{study}/logs/DGE_DTE/kallisto_quant_{read}.log"
        shell:
            "kallisto quant -i {params.kallisto_index} -o {params.outdir} -t {params.threads} {input.fwd} {input.rev} 2> {log}"


else:
    rule kallisto_quant_single:
        input:
            fastq="{study}/reads/{study}/{read}/{read}.fastq.gz",
        output:
            abundance_h5="{study}/DGE_DTE/quant/{read}/abundance.h5",
            abundance_tsv="{study}/DGE_DTE/quant/{read}/abundance.tsv",
            run_info="{study}/DGE_DTE/quant/{read}/run_info.json"
        params:
            outdir="{study}/DGE_DTE/quant/{read}/",
            fragment_length=config["fragment_length"],
            fragment_length_sd=config["fragment_length_sd"],
            kallisto_index=expand("{study}/DGE_DTE/index/kallisto_index.idx", study=studies),
            threads=8
        conda:
            "./envs/kallisto.yaml"
        log:
            "{study}/logs/DGE_DTE/kallisto_quant_{read}.log"
        shell:
            "kallisto quant -i {params.kallisto_index} -o {params.outdir} -t {params.threads} --single {input.fastq} -l {params.fragment_length} -s {params.fragment_length_sd} 2> {log}"

rule star_index:
    output:
        star_index=expand("{study}/Alignment/star_index/SAindex", study=studies)
    priority: 1
    conda:
        "./envs/star.yaml"
    params:
        fasta= config["fasta"],
        gtf_AS_isoforms= config["gtf_AS_isoforms"],
        outdir=expand("{study}/Alignment/star_index/", study=studies)
    log: expand("{study}/logs/Alignment/star_index{study}Log.out", study=studies)
    shell:
        """
        STAR --runThreadN 6 \
        --runMode genomeGenerate \
        --genomeDir {params.outdir} \
        --genomeFastaFiles {params.fasta} \
        --sjdbGTFfile {params.gtf_AS_isoforms} \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 10 \
        2> {log}
        """


rule whippet_sig_events:
    input:
        diff="{study}/AS/DeltaPSI/{comparison_whippet}.diff.gz"
    output:
        diff_unzipped="{study}/AS/DeltaPSI/global/{comparison_whippet}.diff",
        diff_sig_global="{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff"
    params:
        AS_DeltaPSI_threshold=config["AS_Whippet_DeltaPSI_threshold"],
        AS_Whippet_Probability_threshold=config["AS_Whippet_Probability_threshold"]
    shell:
        """
        gzip -dc < {input.diff} > {output.diff_unzipped}
        awk "{{ if ( (\$8 >= {params.AS_DeltaPSI_threshold} || \$8 <= -{params.AS_DeltaPSI_threshold}) && \$9 >= {params.AS_Whippet_Probability_threshold} ) {{print}} }}" {output.diff_unzipped} > {output.diff_sig_global}
        """
#percentage of mapping from whippet .map.gz:

rule whippet_uncompress:
    input:
        whippet_mapping="{study}/AS/{read}.map.gz"
    output:
        whippet_mapping_dc="{study}/AS/{read}.map"
    shell:
        """
        gzip -dc < {input.whippet_mapping} > {output.whippet_mapping_dc}
        """

rule whippet_mapping_summary:
    input:
        whippet_mapping_dc=expand("{study}/AS/{read}.map", study=studies, read=read_files)
    output:
        whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies)
    conda:
        "./envs/whippet_summary.yaml"
    script:
        "./scripts/whippet_mapping_summary.py"



#possibly best to aggregate them in python into one: this could be done adding an extra column for the condition
rule whippet_AS_summary:
    input:
        event_tables=expand("{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff", study=studies, comparison_whippet=comparison_whippet_wc)
    output:
        event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies)
    params:
        tax_id=tax_id
    conda:
        "./envs/whippet_summary.yaml"
    script:
        "./scripts/whippet_AS_global_summary.py"


rule whippet_summary_plots:
    input:
        whippet_mapping_summary=expand("{study}/AS/Summary/Whippet_mapping_summary.csv", study=studies),
        event_table=expand("{study}/AS/Summary/summary_sig_global.csv", study=studies)
    output:
        whippet_mapping_summary_plot=expand("{study}/AS/Summary/Whippet_mapping_summary_plot.svg", study=studies),
        whippet_AS_type_frequency_plot=expand("{study}/AS/Summary/Whippet_Frequency_AS_Event_Types_plot.svg", study=studies),
        whippet_DAS_gene_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_gene_overlap_among_comparisons_plot.svg", study=studies),
        whippet_DAS_event_overlap_among_comparisons_plot=expand("{study}/AS/Summary/whippet_DAS_event_overlap_among_comparisons_plot.svg", study=studies)
    conda:
        "./envs/r_plots.yaml"
    script:
        "./scripts/AS_Whippet_summary.R"
        
rule dge:
    input:
        abundance_h5=expand("{study}/DGE_DTE/quant/{run}/abundance.h5", study=studies, run=read_files) #dge_runs so remove replicates
    output:
        total_heatmap=expand("{study}/DGE_DTE/Summary/total_heatmap.png", study=studies),
        pca=expand("{study}/DGE_DTE/Summary/pca.png", study=studies)
    params:
        dge_names=dge_names,
        gtf_AS_isoforms= config["gtf_AS_isoforms"]
    conda:
        "./envs/deseq2.yaml"
    script:
        "./scripts/deseq2_pairwise.R"



if config["gene_ids_of_interest"]:
    rule DeltaPSI_gene_ids_of_interest:
        input:
            diff_sig_global="{study}/AS/DeltaPSI/global/{comparison_whippet}_sig_global.diff"
        output:
            diff_sig_goi="{study}/AS/DeltaPSI/genes_of_interest/{comparison_whippet}_sig_goi.diff"
        params:
            goi="|".join(list(config["gene_ids_of_interest"].split(",")))
        log:
            "{study}/logs/diff_as_{comparison_whippet}.log"
        shell: #use head -n 1 instead of awk as awk prints first element of first row only
            """
            ( awk "NR==1{{print}}" {input.diff_sig_global}; awk "/{params.goi}/{{print}}" {input.diff_sig_global} ) > {output.diff_sig_goi}
            """

            
    rule bam_file_list:
        input:
            #bam_list=expand("{study}/AS/BAM_sorted/{read}.bam", study=studies, read = read_files)
            input_bam_A=input_bam_A,
            input_bam_B=input_bam_B
            #lambda wildcards : expand("{study}/AS/BAM_sorted/{sample}.bam", sample= comparison_dict[wildcards.comparison_whippet]["A"], study=studies),
            #lambda wildcards : expand("{study}/AS/BAM_sorted/{sample}.bam", sample= comparison_dict[wildcards.comparison_whippet]["B"], study=studies)
            #whippet_out_sorted="{study}/AS/BAM_sorted/{read}.bam"
        params:
            comparison_dict_bam=comparison_dict_v2,
            comparison_name= lambda wildcards : wildcards.comparison_whippet
        conda:
            "./envs/comparison_bam_list.yaml"
        output:
            bam_tsv="{study}/AS/DeltaPSI/{comparison_whippet}_bam.tsv" #names for bam files will also change, maybe it does not need to be included.

        script:
            "./scripts/comparison_bam_list.py"
            
    rule goi_full_coordinates:
        output:
            gtf_genes="{study}/AS/DeltaPSI/genes_of_interest/gtf_genes_coordinates.txt"
        params:
            gtf_AS_isoforms= config["gtf_AS_isoforms"]
        conda:
            "./envs/GenomicFeatures.yaml"
        script:
            "./scripts/max_coords_AS_event_genes.R"

    rule plottable_AS_events:
        input:
            diff_sig_goi="{study}/AS/DeltaPSI/genes_of_interest/{comparison_whippet}_sig_goi.diff",
            bam_tsv="{study}/AS/DeltaPSI/{comparison_whippet}_bam.tsv",
            gtf_genes="{study}/AS/DeltaPSI/genes_of_interest/gtf_genes_coordinates.txt"
        output:
            mock_plot="{study}/AS/sashimi_plots/{comparison_whippet}.txt"
        conda:
            "./envs/ggsashimi.yaml"
        params:
            gtf_AS_isoforms= config["gtf_AS_isoforms"],
            comparison_name= lambda wildcards : wildcards.comparison_whippet,
            study_name=lambda wildcards : wildcards.study
        script:
            "./scripts/AS_goi_plotting.py"
    
    rule plottable_AS_events_detail:
        input:
            diff_sig_goi="{study}/AS/DeltaPSI/genes_of_interest/{comparison_whippet}_sig_goi.diff",
            bam_tsv="{study}/AS/DeltaPSI/{comparison_whippet}_bam.tsv",
            gtf_genes="{study}/AS/DeltaPSI/genes_of_interest/gtf_genes_coordinates.txt"
        output:
            mock_plot_detail="{study}/AS/sashimi_plots/{comparison_whippet}_detail.txt"
        conda:
            "./envs/ggsashimi.yaml"
        params:
            gtf_AS_isoforms= config["gtf_AS_isoforms"],
            comparison_name= lambda wildcards : wildcards.comparison_whippet,
            study_name=lambda wildcards : wildcards.study
        script:
            "./scripts/AS_goi_plotting_detail.py"
            
    rule dge_geneplots:
        input:
            abundance_h5=expand("{study}/DGE_DTE/quant/{run}/abundance.h5", study=studies, run=read_files) #dge_runs to remove replicates
        output:
            mock_gene_plot="{study}/DGE_DTE/geneplots.txt",
            goi_heatmap_all_cond="{study}/DGE_DTE/Summary/goi_all_cond_heatmap.png"
        params:
            goi=config["gene_ids_of_interest"],
            dge_names=dge_names,
            gff_dge= config["gff_dge"],
            study_name=lambda wildcards : wildcards.study
        conda:
            "./envs/deseq2.yaml"
        script:
            "./scripts/deseq2_plotcounts.R"

