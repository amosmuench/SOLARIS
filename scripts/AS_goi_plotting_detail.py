import pandas as pd
import os
from sys import exit
diff_sig_goi=snakemake.input["diff_sig_goi"]
bam_list=snakemake.input["bam_tsv"]
#gtf_genes=snakemake.input["gtf_genes"]
gtf_AS_isoforms=snakemake.params["gtf_AS_isoforms"]
comparison_name=snakemake.params["comparison_name"]
study_name=snakemake.params["study_name"]


diff_goi_df=pd.read_csv(f"{snakemake.input[0]}", header=0, sep="\t", index_col=False )
print("DIff goi table", diff_goi_df)

diff_goi_coord=diff_goi_df["Coord"].tolist()
chrom=[i.split(':', 1)[0] for i in diff_goi_coord]
coords=[i.split(':', 1)[1] for i in diff_goi_coord]
first_coord_expanded=[i.split('-', 1)[0] for i in coords]
first_coord_expanded = [int(i) for i in first_coord_expanded]
first_coord_expanded = [str(i-100) for i in first_coord_expanded] # add basepairs to get flanks


second_coord_expanded=[i.split('-', 1)[1] for i in coords]
second_coord_expanded=[int(i) for i in second_coord_expanded]
second_coord_expanded = [str(i+100) for i in second_coord_expanded]
#diff_goi_coord=list(map(lambda i, j, h : i + ":"+ j + "-" + h, zip(coords, first_coord_expanded, second_coord_expanded) )) 
coords_merged=[i + "-"+ j for i,j in zip(first_coord_expanded, second_coord_expanded)]
diff_goi_coord=[i + ":"+ j for i,j in zip(chrom, coords_merged)]
print("coord changes")
print(diff_goi_coord)
diff_goi_name=diff_goi_df["Gene"].tolist()





#gtf_genes=pd.read_csv(f"{gtf_genes}", header=0, sep=",", index_col=0 )
#gene_id
#Chr3:19864115-19864153

 #seqnames  start    end strand   gene_id
 #  Chr1   3631   5899      + AT1G01010

#gtf_genes['Coord'] = gtf_genes.agg(lambda x: f"{x['seqnames']}:{x['start']}-{x['end']}", axis=1)

#search = ['FR-001', 'FR-002', 'FR-003', 'FR-004']
#search=diff_goi_name
#df['FR'] = df['Description'].str.findall('(' + '|'.join(search) + ')')
#print(gtf_genes)

#gtf_genes['hits'] = gtf_genes['gene_id'].str.findall('(' + '|'.join(diff_goi_name) + ')')
#print(gtf_genes)
#gtf_genes = gtf_genes[gtf_genes['hits'].astype(bool)]
#print(gtf_genes)
#diff_goi_coord=diff_goi_df["Coord"].tolist()
#diff_goi_name=diff_goi_df["Gene"].tolist()
#get entire gene coordinates from table  using findall with diff_goi name as query

try:
    os.makedirs(f"{study_name}" + "/AS/sashimi_plots" + "/"+ "detailed" +"/" + f"{comparison_name}")
except FileExistsError:
    # directory already exists
    pass

outputname_list=[]
for i, j  in zip(diff_goi_name, diff_goi_coord):
    outputname=f"{study_name}" + "/AS/sashimi_plots" + "/" + "detailed" +"/"+ f"{comparison_name}" +"/" + f"{comparison_name}" + "_" + f"{i}" + "_" + f"{j}" + "detailed" 
    outputname_list.append(outputname)
    with open(f"{snakemake.output[0]}", "w") as f:
        f.write(outputname)
if len(diff_goi_df.index)==0:
  print("nothing to plot")
  with open(f"{snakemake.output[0]}", "w") as f:
    f.write("no events found")
  exit(0)
else:

  cwd=os.getcwd()

  path_to_script=cwd + "/scripts/sashimi-plot.py"

  for j in range(0, len(diff_goi_coord)):
    os.system(f"python {path_to_script} -c {diff_goi_coord[j]} -F svg -g {gtf_AS_isoforms} -b {bam_list} -o {outputname_list[j]} --ann-height 5")
     # -M minimum coverage to draw read junction, set to 10 maybe
#best would be to get a dynamic value for the ann height depending on the amount of transcripts for a gene in the gtf
