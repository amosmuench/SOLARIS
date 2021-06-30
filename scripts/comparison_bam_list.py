
import csv
import os
#bam_list=snakemake.input["bam_list"]
a_bam=snakemake.input["input_bam_A"]
b_bam=snakemake.input["input_bam_B"]

print("b-bam")
print(b_bam)


  
   
#print(bam_list)
comparison_dict_bam=snakemake.params["comparison_dict_bam"]
comparison_name=snakemake.params["comparison_name"]
#print(comparison_dict_bam)
#print(comparison_name)
#whippet_out_sorted="{study}/AS/BAM_sorted/{read}.bam"
print("dict")
print(comparison_dict_bam)
a_list=comparison_dict_bam[f"{comparison_name}"]["A"]
b_list=comparison_dict_bam[f"{comparison_name}"]["B"]
 
    #problem with file links. currently, a list of paths to bam files gets imported.
    #get first element, os.path.dirname  assign to directory var this is so dirty I want to wash my hands

path_to_bam=[]
id_bam=[]
#directory=os.path.dirname(a_bam[0])
directory=os.getcwd()#+ "/"+directory
print(directory)
print(a_list)
print(b_list)
a_id, old_a_bam =zip(*[a_list])
a_id=list(a_id)
#a_bam=list(a_bam)
b_id, old_b_bam =zip(*[b_list])
b_id=list(b_id)
#b_bam=list(b_bam)
#if any(isinstance(a_bam, list)):
 # rep=len(a_bam[0])
  #a_id=[i for j in a_id for i in (j,)*rep]
  
  
#if any(isinstance(el, list) for el in a_bam):
  #rep=len(a_bam[0])
  #print("reps", rep)
  
 # a_bam = [item for sublist in a_bam for item in sublist]
  
  #a_list=comparison_dict_bam[f"{comparison_name}"]["A"]
  #a_list=[i for j in a_list for i in (j,)*rep]
  #print(a_list)
  #print(a_bam)
  
#if any(isinstance(el, list) for el in a_bam):
  #rep=len(a_bam[0])
  #print("reps", rep)
  
  #a_bam = [item for sublist in a_bam for item in sublist]
  
  #a_list=comparison_dict_bam[f"{comparison_name}"]["A"]
  #a_list=[i for j in a_list for i in (j,)*rep]
 # print(a_list)
  #print(a_bam)



#if any(isinstance(el, list) for el in b_bam):
  #rep=len(a_bam[0])
  #print("reps", rep)
  
  #b_bam = [item for sublist in b_bam for item in sublist]
  
  #b_list=comparison_dict_bam[f"{comparison_name}"]["B"]
  #b_list=[i for j in b_list for i in (j,)*rep]
  #print(b_list)
  #print(b_bam)
print("a_id")
if any(isinstance(el, list) for el in a_id):
  a_id= [item for sublist in a_id for item in sublist]
print(a_id)
if any(isinstance(el, list) for el in b_id):
  b_id= [item for sublist in b_id for item in sublist]
print(b_id)
for i in range(0, len(a_bam)): #this will break if multiple in list
    print("a")
    id_bam.append(a_id[i])
    bam_a = os.path.join(directory, f"{a_bam[i]}")
    path_to_bam.append(bam_a)
for j in range(0, len(b_bam)):
    print("b")
    print(b_id[j])
    id_bam.append(b_id[j])
    bam_b = os.path.join(directory, f"{b_bam[j]}")
    print(bam_b)
    path_to_bam.append(bam_b)


print(id_bam)

print(snakemake.output["bam_tsv"])
with open(f"{snakemake.output[0]}", "w") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(zip(id_bam, path_to_bam))

    #bam_tsv=dict(zip(id_bam, path_to_bam))
    #print(bam_tsv)
    #bam_tsv=pd.Series(bam_tsv).to_frame()
    #bam_tsv=pd.DataFrame.from_records(bam_tsv)
    #bam_tsv.to_csv(snakemake.output["bam_tsv"], index=False, header=False, sep="\t")
