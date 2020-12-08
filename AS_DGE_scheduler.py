#modules
import requests
import os
#import xmltodict
from io import StringIO
import pandas as pd
pd.options.mode.chained_assignment = None #to supress warning, as far as I know the chained_assignment is not applicable to the script
from urllib.request import urlopen
import xml.etree.ElementTree as et
from itertools import combinations, product
from itertools import chain
import itertools
import sys
import pickle


#best would be a class:
#class RNA_Seq_Study:
    #class attribute
    #sequencing_object="RNA-Seq Study"
    #def __init__(self, study_name, seq_type, read_length, sample_accession, run_list):
    #    self.study_name=study_name
    #    self.type=seq_type
    #    self.read_length=read_length
    ##    self.sample_accession=sample_accession
    #    self.run_list=run_list
    #    #add contrasts
    #def __str__(self):
    #    return f"{self.study_name} is a {self.read_length} bp {self.type} {RNA_Seq_Study.sequencing_object} with {len(self.run_list)} Runs."

    #def ...
    #add functions as member functions...


#pickle instance of the class, import class, load into snakefile


#functions

def is_redundant(df, A, B):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    left=A
    right=B
    if len(df.groupby(f'{left}')) == len(df.groupby([f'{left}', f'{right}'])):
        return True
    elif len(df.groupby(f'{right}')) == len(df.groupby([f'{right}', f'{left}'])):
        return True
    else:
        return False


def redundant_groups_finder(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]

    cols = [c for c in df.columns]

    redundant_groups = []
    if len(cols)==1:
        return(redundant_groups)
    else:
        index_l = 0
        while index_l < len(cols)-1:
            group_n = []
            index_r = index_l+1
            while index_r < len(cols):
                if is_redundant(df, cols[index_l], cols[index_r]):
                #print( cols[index_l], cols[idx_right])
                    group_n.append(cols.pop(index_r))
                else:
                    index_r += 1
            if group_n:
                redundant_groups.append(group_n + [cols[index_l]])
            index_l += 1

            return(redundant_groups)

def read_type_det(read_type): #best would be one var to be returned
    for i in read_type:
        if len(i) ==2:
            paired = True
            single = False
            #print("paired")
        else:
            #print("single")
            paired = False
            single = True
    return(single, paired)



def factor_levels_dict(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    factor_list=df.columns.values.tolist()
    factors_levels= {col: list(df[col].unique()) for col in factor_list}

    return(factors_levels)


def factors_levels_size(dict):

    all_levels=list(dict.values())
    level_size=[]
    for i in all_levels:
        level_size.append(len(i))
    return(level_size)



def multiplier(list):
    multiplication_result=1
    for i in list:
        temp=multiplication_result*i
        multiplication_result=temp
    return(multiplication_result)

#this is an ugly argument
def replicates_def(df, level_size):
    full_df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    level_size=level_size
    replicates=None
    number_of_replicates=0
    if len(full_df.index) == multiplier(level_size):

        replicates=False
        number_of_replicates=1

        return(replicates, number_of_replicates)

    else:
        number_of_replicates=int(len(full_df.index)/multiplier(level_size))

        replicates=True

    return(replicates, number_of_replicates)


def pairs(*lists):
    for i in combinations(lists, 2):
        for pair in product(*t):
            if pair[0] != pair[1]:
                yield pair

def unique_everseen_custom(it):
    checked = set()
    checked_added = seen.add
    for tuple in it:
        if (tuple not in checked) and (tuple[::-1] not in checked):
            checked_added(tuple)
            checked_added(tuple[::-1])
            yield tuple

def get_constrast_indices(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    length=len(df.index)
    a = range(0,length,1)
    #b = range(0,len(df.index),1)
    list_ = set(pairs(a, a))
    res = list(unique_everseen_custom(list_))
    contrasts_indices = sorted(res, key=lambda tup: tup[0])
    print(contrasts_indices)
    return(contrasts_indices)

def factor_with_highest_level(factors_levels):
    level_length=[]
    length_of_highest_level=None
    for i in list(factors_levels.values()):
        level_length.append(len(i))
    length_of_highest_level=max(level_length)
    index=-1
    for i in list(factors_levels.values()):
        index+=1
        if len(i) == max(level_length):
            break
        else:
            continue
    highest_level_factor=list(factors_levels.keys())[index]
    return(highest_level_factor)


def subset_sample_table(subset_by_level, df):

    subset_column=None
    subsetted_df={}
    dict=factor_levels_dict(df)
    subset_column=[factor for factor in dict.keys() if f'{subset_by_level}' in dict[factor]]
    subset_column="".join(subset_column)
    subsetted_df=df[df[f"{subset_column}"]==subset_by_level]
    #subsetted_df=df.query(f'{subset_column}=={subset_by_level}')
    return(subsetted_df)


def replicates_merger(df, non_unique_col_list, replicate_nr):
    keep_list=['sample_accession', 'run_accession', 'read_file']
    a=non_unique_col_list[0] #get first element, this is questionable, maybe the first element will not be the unique one
    df_temp=(df.groupby(df[a].tolist())
           .apply(lambda x: tuple(x.index))
           .reset_index(name='repl_index'))

    repl_index=df_temp['repl_index'].tolist()
    list_index=[]
    for i in repl_index:
        list_index.append(list(i))

    sample_accession_new=[]
    run_accession_new=[]
    read_file_new=[]
    for i in list_index:
        sample_accession_new_row=[]
        sample_accession_new_row.append(df['sample_accession'].iloc[i].tolist())
        sample_accession_new_row=list(chain.from_iterable(sample_accession_new_row))
        sample_accession_new.append(sample_accession_new_row)

        run_accession_new_row=[]
        run_accession_new_row.append(df['run_accession'].iloc[i].tolist())
        run_accession_new_row=list(chain.from_iterable(run_accession_new_row))
        run_accession_new.append(run_accession_new_row)

        read_file_new_row=[]
        read_file_new_row.append(df['read_file'].iloc[i].tolist())
        read_file_new_row=list(chain.from_iterable(read_file_new_row))
        read_file_new.append(read_file_new_row)


    replicates_list=list(range(0, replicate_nr))
    replicates_columns=[]
    for i in replicates_list:
        replicates_columns.append(str(i))
    naming=["sample_accession_repl", "run_accession_repl","read_file_repl"] #3 three times
    out=list(itertools.product(naming, replicates_columns))
    detouple=()
    for i in out:
        detouple+=i
    detouple=iter(detouple)
    replicates_columns=[i+next(detouple, '') for i in detouple]

    sample_accession_new=pd.DataFrame(sample_accession_new)
    run_accession_new=pd.DataFrame(run_accession_new)

    read_file_new=pd.DataFrame(read_file_new)

    out_df=pd.concat([sample_accession_new, run_accession_new, read_file_new], axis=1)

    out_df.columns=replicates_columns
    out_df["sample_accession"]=out_df.filter(regex=("sample_accession*")).values.tolist()
    out_df["run_accession"]=out_df.filter(regex=("run_accession*")).values.tolist()
    out_df["read_file"]=out_df.filter(regex=("read_file*")).values.tolist()
    #out_df.filter(regex=("d.*"))
    out_df.drop(replicates_columns, axis=1, inplace=True)
    df=df.drop(keep_list, axis=1)
    df=df.drop_duplicates(keep='first', inplace=False)
    df=df.reset_index(drop=True)
    out_df=pd.concat([df, out_df], axis=1)
    return(out_df)

#function runs recursively, this is probably quite slow, better: dp approachh
def drop_redundant_v2(df, redundant_groups):
    redundant_groups = [item for sublist in redundant_groups for item in sublist]
    redundant_groups.sort()
    if not redundant_groups:
        return(df)
    elif len(df.columns)==1:
        return(df)
    else:
        df = df.reindex(sorted(df.columns), axis=1)
        redundant_groups=df.columns.tolist()
        to_be_merged_df=df
        red_pairs=itertools.combinations(redundant_groups, r=2)
        red_pairs_list=[]
        for i in red_pairs: #unnecessary
            temp=tuple(i)
            red_pairs_list.append(temp)
        for i in red_pairs_list:
            int1=to_be_merged_df.groupby(f'{i[0]}')
            if len(to_be_merged_df.columns)>2:
                int2=to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}'])
                grpby_first_sec_count=to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}']).count().iloc[:,0].tolist()[0]
            else:
                grpby_first_sec_count=to_be_merged_df.groupby([f'{i[1]}']).count().iloc[:,0].tolist()[0]
            grpby_first_count=to_be_merged_df.groupby(f'{i[0]}').count().iloc[:,0].tolist()[0]

            grpby_first_unique=to_be_merged_df.groupby(f'{i[0]}').nunique().iloc[:,0].tolist()[0]



            if int(grpby_first_count/grpby_first_unique)==grpby_first_sec_count:

                to_be_merged_df.drop([f'{i[0]}'], axis=1, inplace=True)

                new_red_group_member=''.join(i)

                new_redundant_groups=redundant_groups
                new_redundant_groups.remove(i[0])
                new_redundant_groups.remove(i[1])
                new_redundant_groups.append(new_red_group_member)

                return(drop_redundant_v2(to_be_merged_df, new_redundant_groups))
                continue
                #return
                #return(to_be_merged_df)
                #this needs to run again with only grouping by remaining vs remaining + merged

            elif len(to_be_merged_df.groupby(f'{i[0]}')) == len(to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}'])):
                 to_be_merged_df.drop([f'{i[0]}'], axis=1, inplace=True)
                 return(df)
            else:
                return(df)

def detect_multiple_studies(df):
    for y, x in df.groupby(df.isnull().dot(df.columns)):

        substudies = {x : y.dropna(1) for x, y in df.groupby(df.isnull().dot(df.columns))}

        #this is still inversed!!
        #how to reverse??
    print("Substudies found:")
    substudies_names=list(reversed(substudies.keys()))
    print(substudies_names)
    for k, v in substudies.items():
        print("Substudy that does not contain: ", k, "\n", v )
    substudy=input("Choose name of substudy to excluded from the analysis: ")
    if substudy in substudies:
        substudy=dict(substudies[substudy])
        df=pd.DataFrame(substudy)
        return(df)

def contrast_prompter():

    global sample_attributes
    if sample_attributes.isnull().values.any():
        sample_attributes=detect_multiple_studies(sample_attributes)
        sample_attributes.reset_index(drop=True, inplace=True)

    else:
        pass
    redundant_groups=redundant_groups_finder(sample_attributes)
    if redundant_groups:
        [redundant_groups] = redundant_groups
        redundant_df=sample_attributes[redundant_groups]
        redundant_df_merged=drop_redundant_v2(redundant_df,redundant_groups)
        rest=[col for col in sample_attributes.columns if col not in redundant_groups]
        rest=sample_attributes[rest]
        sample_attributes=pd.concat([redundant_df_merged, rest], axis=1)

        factors_levels=factor_levels_dict(sample_attributes)
        level_size=factors_levels_size(factors_levels)
        replicates_nr=replicates_def(study_table, level_size)
    else:
        factors_levels=factor_levels_dict(sample_attributes)
        level_size=factors_levels_size(factors_levels)
        replicates_nr=replicates_def(study_table, level_size)
    if replicates_nr[0]: #true there are replicates
        sample_attributes=replicates_merger(sample_attributes, list(factors_levels.keys()), replicates_nr[1])
    else:
        sample_attributes=sample_attributes


    print(f"{study} ")
    read_type_det
    global paired #no good practice to set as global arguemnt missing somewhere in function
    global single
    if paired:
        print("is paired-end sequencing data")
    else:
        print("is single-end sequencing data")
    if replicates_nr[1]==0:
        print(f"with {len(sample_attributes.index)} samples")
    else:
        print(f"with {len(sample_attributes.index)*replicates_nr[1]} samples") #this needs to be length of sample attributes before replicates merged: len df * rep nr
    if replicates_nr[0]:
        print(f"and contains replicates.")
    else:
        print("without replicates.")
    #print(f"and contains replicates.")
    print(f"Number of Replicates: {replicates_nr[1]}.")
    print(f"Experimental factors are {list(factors_levels.keys())}")
    print(f"with levels {list(factors_levels.values())}")
    print("\nSample table:")

    print(sample_attributes[sample_attributes.columns[~sample_attributes.columns.isin(sample_info_col_list)]])
    print(f"\nPipeline will investigate Differential Gene/Transcript Expression (DGE/DTE) and \nnovel and annotated Alternative Splicing (AS) events\nusing pairwise comparisons between all conditions.\n")
#missing: load from yaml file which gff and  fasta will be used
    print("Default settings will be used: \nGenome Annotation for DGE/DTE: AtRTD2.gff \nGenome Annotation for AS: AtRTD2_QUASI.gff \nGenomic Sequence: TAIR10.fasta \nDESeq2 using factors to build GLM, Wald-Test for DGE/DTE. \nWhippet detection of AS events DeltaPSI >=0.1, Probability >=0.9")
    contrasts_indices=get_constrast_indices(sample_attributes)

    not_started=True
    while not_started:
        if len(contrasts_indices) >=25:
            reduction_needed=True
        else:
            reduction_needed=False
        while reduction_needed:
            print(f"The number of pairwise comparisions exceeds 25 (is: {len(contrasts_indices)}).")
            print(f"Hence, a factor needs to be removed.")
            highest_level_factor=factor_with_highest_level(factors_levels)
            print(f"The factor with the highest amount of levels is: {highest_level_factor}")
            subset_by_level=input("Choose level of factor by which the Study should be subsetted: ")
        #check whether is in
            confirm_subset_by_level=input(f"Study will be subsetted by level {subset_by_level} enter y to continue, n to abort: ")
            if confirm_subset_by_level=="y":
            #sample_attributes_temp=sample_attributes.copy()
                subsetted_sample_attributes=subset_sample_table(subset_by_level, sample_attributes)
            #this raises a key error if sample_attributes is not the none_reduant version

                reduced_contrasts_indices=get_constrast_indices(subsetted_sample_attributes)

                print(f"new number of pairwise comparisons: {len(reduced_contrasts_indices)}")
                if len(reduced_contrasts_indices)<25:
                    reduction_needed=False
                    sample_attributes=subsetted_sample_attributes
                    contrasts_indices=reduced_contrasts_indices



                else:
                    print("contrast number still to high")
                    sample_attributes=subsetted_sample_attributes
                    contrasts_indices=reduced_contrasts_indices



            else:
                return

        print("Sample Table which will be used for pairwise comparisons:")
        print(sample_attributes)
        print("\n")
        #print("Contrasts indices:")
        sample_attributes_confirmed=sample_attributes
        contrast_indices_confirmed=contrasts_indices
        sample_attributes_confirmed.reset_index(drop=True, inplace=True)
        comp_A=[]
        comp_B=[]
        #for i in contrast_indices_confirmed:
        #    comp_A.append(sample_attributes_confirmed.iloc[i[0]])
        #    comp_B.append(sample_attributes_confirmed.iloc[i[1]])

        #df[sample_attributes_confirmed.columns[~sample_attributes_confirmed.columns.isin(sample_info_col_list)]]

        sample_attributes_confirmed['comp'] = sample_attributes_confirmed[sample_attributes_confirmed.columns[~sample_attributes_confirmed.columns.isin(sample_info_col_list)]].apply(
            lambda x: '_'.join(x.astype(str)),
            axis=1
        )
        comp=sample_attributes_confirmed['comp'].tolist()
        run_list=sample_attributes_confirmed['read_file'].tolist()

        run_accession_list=sample_attributes_confirmed['run_accession'].tolist()

        sample_acc_list=sample_attributes_confirmed['sample_accession'].tolist()

        A=[i[0] for i in contrast_indices_confirmed]
        comp_A_names=[comp[i] for i in A]
        print("First pair of comparison: ", comp_A_names)
        print("\n")
        run_list_A=[run_list[i] for i in A]
        B=[i[1] for i in contrast_indices_confirmed]
        comp_B_names=[comp[i] for i in B]

        print("Second pair of comparison:", comp_B_names)
        run_list_B=[run_list[i] for i in B]

        sample_acc_list_A=[sample_acc_list[i] for i in A]
        sample_acc_list_B=[sample_acc_list[i] for i in B]


        run_accession_list_A=[run_accession_list[i] for i in A]
        run_accession_list_B=[run_accession_list[i] for i in B]
        ################ names ############
        contrast_names = [i + "_vs_"+ j for i, j in zip(comp_A_names, comp_B_names)]
        contrast_names=[i.replace(" ", "") for i in contrast_names]

        print(contrast_names)

        print(replicates_nr[1])
        print(range(0,replicates_nr[1]))
        replicate_suffix=list(range(0,replicates_nr[1]))*len(comp_A_names)
        comp_A_names_replicates=[i+ "_R" + str(j) for i,j in zip(comp_A_names,replicate_suffix)]
        cond_A_acc=zip(comp_A_names_replicates, run_accession_list_A)


        comp_B_names_replicates=[i+ "_R" + str(j) for i,j in zip(comp_B_names,replicate_suffix)]
        print(comp_A_names_replicates)
        run_accession_list_B
        cond_B_acc=zip(comp_B_names_replicates, run_accession_list_B)

        key_1 = 'A'
        key_2 = 'B'
        comparison_dict_v2={k : {key_1 : v1, key_2 : v2} for k,v1,v2 in zip(contrast_names, cond_A_acc, cond_B_acc)}

        A_list=list(zip(comp_A_names, sample_acc_list_A, run_list_A))
        #print("A_list", A_list)
        B_list=list(zip(comp_B_names, sample_acc_list_B, run_list_B))
        complete_list=list(zip(A_list, B_list))

        dic_A={}
        index=0
        for i in comp_A_names:
            dic_A[i]=run_list_A[index]
            index+=1

        dic_B={}
        index=0
        for i in comp_B_names:
            dic_B[i]=run_list_B[index]
            index+=1

        key_1 = 'A'
        key_2 = 'B'
        comparison_dict={k : {key_1 : v1, key_2 : v2} for k,v1,v2 in zip(contrast_names, run_accession_list_A, run_accession_list_B)}

        print("Number of comparisons:", len(contrasts_indices))
        print("\n")
        print(replicates_nr[0])
        start_pipeline=input("To issue pipleline using pairwise comparison between each to HPC enter y, to abort, enter n: ")
        #when this becomes a class just one instance needs to be pickles
        #some of them are
        if start_pipeline=="y":
            study_total=open(output_dir +'study_total.pkl', 'wb')
            replicates=replicates_nr[0]
            pickle.dump(study, study_total)
            pickle.dump(run_accession_list, study_total)
            pickle.dump(paired, study_total)
            pickle.dump(replicates, study_total)
            pickle.dump(complete_list, study_total)
            pickle.dump(comparison_dict, study_total)
            pickle.dump(comparison_dict_v2, study_total)
            pickle.dump(comp_A_names,study_total)
            pickle.dump(sample_acc_list_A, study_total)
            pickle.dump(run_accession_list_A, study_total)
            pickle.dump(run_list_A, study_total)
            pickle.dump(comp_B_names,study_total)
            pickle.dump(sample_acc_list_B, study_total)
            pickle.dump(run_accession_list_B, study_total)
            pickle.dump(run_list_B, study_total)
            study_total.close()

            os.system("sbatch snake.sh")
            print("Pipeline started.")
            return(sample_attributes, contrasts_indices) #unncessary could be void
            #this needs to be delayed to make sure that the subsetting was done or do it in the function
            #works


            not_started=False
        if start_pipeline=="n":
            exit()

study=str(sys.argv[1])

cwd = os.getcwd()
output_dir=cwd+"/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Retrieving information for study {study} from ENA")

print("""\

                                       ._ o o
                                       \_`-)|_                  GIRAFFE PICTURE
                                    ,""       \               AS WE ARE WAITING FOR THE SERVER...
                                  ,"  ## |   ಠ ಠ.
                                ," ##   ,-\__    `.
                              ,"       /     `--._;)
                            ,"     ## /
                          ,"   ##    /


                    """)


#this could be wrapped into function
study_table_request = requests.get(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={study}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true")

#not catching empty/wrong study numbers yet
study_table=study_table_request.content

study_table=str(study_table, 'utf-8')

study_table=StringIO(study_table)
study_table=pd.read_csv(study_table, sep='\t', lineterminator='\n')

sample_accession=study_table['sample_accession'].tolist()
read_files=study_table['sample_accession'].tolist()

fastq_ftp=study_table['fastq_ftp'].tolist()

temp=[]
for i in fastq_ftp:
    temp.append(i.split(";"))





single, paired=read_type_det(temp)


read_type=[]
first=None
last=None
if paired:
    for i in temp:
        first=i[0]
        first=os.path.basename(first)
        last=i[1]
        last=os.path.basename(last)

        both=(first, last)
        read_type.append(both)

else:

    for i in temp:
    #    print(i)
        first=i[0]
        first=os.path.basename(first)
        read_type.append(first)


run_accession=study_table['run_accession'].tolist()

xml_links=[]
for i in sample_accession:
    xml_links.append(f"https://www.ebi.ac.uk/ena/browser/api/xml/{i}")



xml_list= [urlopen(link).read() for link in xml_links]

#xml_file= open("SAMN10606778.xml", "r")
xtree=[et.ElementTree(et.fromstring(xml)) for xml in xml_list]


xroot=[xtree_i.getroot() for xtree_i in xtree]

#missing: taxid pull and pickle
sample_attributes=[]
for i in xroot:
    xheader_inner=[]
    xvalues_inner=[]
    for j in i.findall(".//SAMPLE_ATTRIBUTE/TAG"):
                    #print(i.text)
        xheader_inner.append(j.text)
    for j in i.findall(".//SAMPLE_ATTRIBUTE/VALUE"):
        #print(i.text)
        xvalues_inner.append(j.text)
    attributes_inner=dict(zip(xheader_inner, xvalues_inner))
    sample_attributes.append(attributes_inner)


#this could be done in one run
index=0
for i in sample_attributes:
    i["sample_accession"]=sample_accession[index]
    index+=1

index=0
for i in sample_attributes:
    i["run_accession"]=run_accession[index]
    index+=1
index=0
for i in sample_attributes:
    i["read_file"]=read_type[index]
    index+=1

sample_attributes=pd.DataFrame(sample_attributes)


sample_attributes = sample_attributes[sample_attributes.columns.drop(list(sample_attributes.filter(regex='ENA'))) ]
print(sample_attributes)

sample_info_col_list=['sample_accession', 'run_accession', 'read_file']

columns_useable=sample_attributes.columns[~sample_attributes.columns.isin(sample_info_col_list)]

for col in columns_useable:
    if len(sample_attributes[col].unique()) == 1:
        sample_attributes.drop(col,inplace=True,axis=1)


sample_attributes_confirmed, contrast_indices_confirmed=contrast_prompter() #could be void
