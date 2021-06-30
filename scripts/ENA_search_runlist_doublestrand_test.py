import os
import requests
from io import StringIO
import pandas as pd
import pickle

#paste snakemake study for accession

studies=snakemake.config["studies"]
directory=snakemake.config["directory"]
#studies="PRJNA510535"

url ="https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+ studies +"&result=read_run"
#"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA510535&result=read_run"

html = requests.get(url).content
print(type(html))
print(html)
s=str(html, 'utf-8')

data=StringIO(s)
df=pd.read_csv(data, sep='\t', lineterminator='\n')
print(df)

run_list=df['run_accession'].tolist()



#integrate study path from  snakemake
output_dir=directory+"/"+studies+"/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_file="run_list.p"
pickle_file=open(output_dir + "run_list.p", "wb")
pickle.dump(run_list, pickle_file)

fastq_ftp=df['fastq_ftp'].tolist()

read_type=[]
for i in fastq_ftp:
    read_type.append(i.split(";"))



def read_type_det(read_type):
    for i in read_type:
        if len(i) ==2:
            paired = True
            single = False
        else:
            paired = False
            single = True
    return(single, paired)

read_type_bool=read_type_det(read_type)

#integrate study path from  snakemake
output_dir=directory+"/"+studies+"/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

pickle_file_run_type=open(output_dir + "run_type.p", "wb")
pickle.dump(read_type_bool, pickle_file_run_type)

#add an md5 chekcsum downstream of the pipeline to see if files were downloaded properly
