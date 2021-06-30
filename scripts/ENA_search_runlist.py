#python -m pip install requests
#python -m pip install pandas
#add both to #!/usr/bin/env python
#beautifulsoup4
#html5lib
#lxml
#io

import os
import requests
from io import StringIO
import pandas as pd
import pickle

#paste snakemake study for accession

studies=snakemake.config["studies"]
directory=snakemake.config["directory"]
#studies="SRP156618"
url ="https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+ studies +"&result=read_run"
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