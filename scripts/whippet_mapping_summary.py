import pandas as pd
import os
import re
#add one column with the name of the df aka condition
print(snakemake.input["whippet_mapping_dc"])

def retrieve_mapping_stats(list_of_paths, key1, key2, key3):
    mapping_summary={}
    sample_name_list=[]
    Mapped_Percent_list=[]
    Multimap_Percent_list=[]
    Novel_Junc_Percent_list=[]
    for i in list_of_paths:
        #read in as csv
        print(i)
        sample_name=os.path.basename(i)
        print(sample_name)
        sample_name=re.sub("\.map$", "", sample_name)
        print(sample_name)
        sample_name_list.append(sample_name)
        #f"{sample_name}"={}
        with open(f"./{i}", 'r') as file:
            text = file.read()
        
        text=re.sub(r'\n', ' ', text)    
        text=re.sub(r'\t', ' ', text)    
        hit_key1=re.search(rf'{key1} (\S+)', text)
        hit_key1=hit_key1.group(1)
        hit_key2=re.search(rf'{key2} (\S+)', text)
        hit_key2=hit_key2.group(1)
        hit_key3=re.search(rf'{key3} (\S+)', text)
        hit_key3=hit_key3.group(1)
        Mapped_Percent_list.append(hit_key1)
        Multimap_Percent_list.append(hit_key2)
        Novel_Junc_Percent_list.append(hit_key3)
        
        #mapping_summary[sample_name]={}
        #mapping_summary[sample_name]={"Multimap_Percent":after_key2}
        #f"{sample_name}"[f"key1"]=after_key1
        #f"{sample_name}"[f"key2"]=after_key2
        #f"{sample_name}"[f"key3"]=after_key3
        #mapping_summary.append(f"{sample_name}") #this needs to be a dictionary
        print(mapping_summary)
    mapping_summary={"Run": sample_name_list, "Mapped_Percent":Mapped_Percent_list, "Multimap_Percent":Multimap_Percent_list, "Novel_Junc_Percent":Novel_Junc_Percent_list}    
    mapping_summary=pd.DataFrame.from_records(mapping_summary)    
    mapping_summary.to_csv(snakemake.output[0], index=False)


retrieve_mapping_stats(list(snakemake.input["whippet_mapping_dc"]), "Mapped_Percent", "Multimap_Percent", "Novel_Junc_Percent")
