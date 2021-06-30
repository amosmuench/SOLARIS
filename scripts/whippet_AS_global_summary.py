import pandas as pd
import os
#add one column with the name of the df aka condition

names=["Gene","Node","Coord","Strand","Type","Psi_A","Psi_B","DeltaPsi","Probability","Complexity","Entropy"]

def df_name2column(list_of_dfs):
    print(list_of_dfs)
    for i in list_of_dfs:
        print(i)
        #read in as csv
        df_name=os.path.basename(i)
        #df_name.split(".")[0]
        df_name = df_name.replace(".diff", "")
        print(df_name)


        df=pd.read_csv(f"{i}", header=0, sep="\t", index_col=False )
        #split away / til .diff
        print(df.head(5))
        df["comparison"]=f"{df_name}"
        print(df.head(5))
        yield df

summary_list=df_name2column(list(snakemake.input["event_tables"])) # this will be a list of file paths
print(summary_list)
summary=pd.concat(summary_list)


summary.to_csv(rf'{snakemake.output["event_table"]}', index=False)
