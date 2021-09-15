import requests
import math
import functools
import os
import argparse
from io import StringIO
import pandas as pd
pd.options.mode.chained_assignment = None
from urllib.request import urlopen
import xml.etree.ElementTree as et

from itertools import combinations, product
from itertools import chain
import itertools

import sys
import pickle
parser = argparse.ArgumentParser()
parser.add_argument("-ena",  help="Study Identifier from ENA, e.g. PRJNA176940", nargs=1)
parser.add_argument('-local', help="Path to CSV file with one column for sample name (replicates should be given the same name) and one column containing the path to the read file (pair-end read files should be in the same cell, separated by a comma)",  nargs=1)
args = parser.parse_args()


def read_type_det(read_type):
    for i in read_type:
        if len(i) ==2:
            paired = True
            single = False
        else:
            paired = False
            single = True
    return(single, paired)

def is_redundant(df, A, B):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]

    left=A
    right=B
    #print(df)
    if len(df.groupby(f'{left}')) == len(df.groupby([f'{left}', f'{right}'])):

        return True
    elif len(df.groupby(f'{right}')) == len(df.groupby([f'{right}', f'{left}'])):
        return True
    else:

        return False

def drop_redundant(df, redundant_groups):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]

    redundant_groups=redundant_groups

    redundant_groups.sort()

    if not redundant_groups:

        return(df)
    else:


        df = df.reindex(sorted(df.columns), axis=1)

        red_pairs=itertools.combinations(redundant_groups, r=2)



        red_pairs_list=[]
        for i in red_pairs:
            temp=tuple(i)
            red_pairs_list.append(temp)

        for i in red_pairs_list:

            int1=to_be_merged_df.groupby(f'{i[0]}')

            int2=to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}'])

            grpby_first_count=to_be_merged_df.groupby(f'{i[0]}').count().iloc[:,0].tolist()[0]

            grpby_first_unique=to_be_merged_df.groupby(f'{i[0]}').nunique().iloc[:,0].tolist()[0]
            grpby_first_sec_count=to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}']).count().iloc[:,0].tolist()[0]


            if int(grpby_first_count/grpby_first_unique)==grpby_first_sec_count:

                to_be_merged_df[f'{i[0]}' + f'{i[1]}'] = to_be_merged_df[[f'{i[0]}', f'{i[1]}']].sum(axis=1)

                to_be_merged_df.drop([f'{i[0]}', f'{i[1]}'], axis=1, inplace=True)

                continue


            elif len(to_be_merged_df.groupby(f'{i[0]}')) == len(to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}'])):
                 to_be_merged_df[f'{i[0]}' + f'{i[1]}'] = to_be_merged_df[[f'{i[0]}', f'{i[1]}']].sum(axis=1)

                 to_be_merged_df.drop([f'{i[0]}', f'{i[1]}'], axis=1, inplace=True)


                 return(df)

            else:
                return(df)


#is_redundant(df,'B','D')
#print(df)


def redundant_groups_finder(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    #print(df)
    cols = [c for c in df.columns]

    redundant_groups = []
    if len(cols)==1:

        return(redundant_groups)
    else:
        idx_left = 0
        while idx_left < len(cols)-1:
            new_group = []
            idx_right = idx_left+1
            while idx_right < len(cols):
                if is_redundant(df, cols[idx_left], cols[idx_right]):

                    new_group.append(cols.pop(idx_right))

                else:

                    idx_right += 1
            if new_group:
                redundant_groups.append(new_group + [cols[idx_left]])


            idx_left += 1


            return(redundant_groups)




def factor_levels_dict(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]

    factor_list=df.columns.values.tolist()

    factors_levels= {col: list(df[col].unique()) for col in factor_list}

    return(factors_levels)





def factors_levels_size(dict):
    #takes dict of factor_levels
    all_levels=list(dict.values())
    level_size=[]
    for i in all_levels:
        level_size.append(len(i))


    return(level_size)



def multiplier(list):
    min_factor=min(list)
    print(min_factor)
    print(len(list))
    print("range", range(1 ,len(list)+1))
    for i in range(1 ,len(list)+1):
        print("i",i)
        temp=min_factor*i
        min_factor=temp
    temp=1
    for i in list:
      print(i)
      temp=i*temp


    print("multiplier")
    print(temp)
    return(temp)

#this is an ugly argument
def replicates_def(df, level_size):
    full_df=df[df.columns[~df.columns.isin(sample_info_col_list)]]

    all_levels=list(full_df.keys())


    grouped=df.groupby(all_levels).size()
    number_of_replicates=min(grouped.values.tolist())



    if number_of_replicates > 1:
        replicates=True
    else:
        replicates=False

    return(replicates, number_of_replicates)

#get pairwise combinations by indexes of the full def, this needs to be wrapped into an entire function


def pairs(*lists):
    for t in combinations(lists, 2):
        for pair in product(*t):
            #Don't output pairs containing duplicated elements
            if pair[0] != pair[1]:
                yield pair

def unique_everseen_custom(it):
    seen = set()
    seen_add = seen.add
    for tup in it:
        if (tup not in seen) and (tup[::-1] not in seen):
            seen_add(tup)
            seen_add(tup[::-1])
            yield tup




def get_constrast_indices(df):
    df=df[df.columns[~df.columns.isin(sample_info_col_list)]]
    length=len(df.index)
    a = range(0,length,1)

    list_ = set(pairs(a, a))
    res = list(unique_everseen_custom(list_))

    contrasts_indices = sorted(res, key=lambda tup: tup[0])

    return(contrasts_indices)
    #add name to it?



#get the factor with the highest amount of levels
def factor_with_highest_level(factors_levels):
    level_length=[]
    length_of_highest_level=None
    for i in list(factors_levels.values()):
        level_length.append(len(i))
    length_of_highest_level=max(level_length)
    index=-1
    for i in list(factors_levels.values()):
        #print(i)
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
    print("reo", replicate_nr)
    keep_list=['sample_accession', 'run_accession', 'read_file']
    a=non_unique_col_list[0]

    test=df.groupby(df[a])
    print(test.describe())
    df_temp=(df.groupby(df[a].tolist())
           .apply(lambda x: tuple(x.index))
           .reset_index(name='repl_index'))


    b=non_unique_col_list[len(non_unique_col_list)-1]

    df_temp2=(df.groupby(df[b].tolist())
           .apply(lambda x: tuple(x.index))
           .reset_index(name='repl_index'))


    #print(a)
    repl_index=df_temp['repl_index'].tolist()
    repl_index2=df_temp2['repl_index'].tolist()
    min_repl1=[]
    min_repl2=[]

    for i in range(0,len(repl_index)):
        min_repl1.append(len(repl_index[i]))
    for j in range(0, len(repl_index2)):
        min_repl2.append(len(repl_index2[j]))
    indices1 = [i for i, x in enumerate(min_repl1) if x == min(min_repl1)]
    max1 = [i for i, x in enumerate(min_repl1) if x == max(min_repl1)]
    indices2=[i for i, x in enumerate(min_repl2) if x == min(min_repl2)]
    max2 = [i for i, x in enumerate(min_repl2) if x == max(min_repl2)]

    max1=[repl_index[i] for i in max1]
    max2=[repl_index2[i] for i in max2]

    max1=[item for t in max1 for item in t]
    max2=[item for t in max2 for item in t]
    print(max1)
    print(max2)
    add=list(set(max1) & set(list(max2) ))
    add=tuple(add)
    print(add)
    repl_index=[repl_index[i] for i in indices1]
    print("r1")
    print(repl_index)
    repl_index2=[repl_index2[i] for i in indices2]
    print("r")
    print(repl_index2)


    repl_index=repl_index#+repl_index2+[add]

    print(repl_index)
    list_index=[]
    for i in repl_index:
        list_index.append(list(i))


    sample_accession_new=[]
    run_accession_new=[]
    read_file_new=[]
    print(list_index)
    for i in list_index:
        sample_accession_new_row=[]
        sample_accession_new_row.append(df['sample_accession'].iloc[i].tolist())
        sample_accession_new_row=list(chain.from_iterable(sample_accession_new_row))
        sample_accession_new.append(sample_accession_new_row)

        run_accession_new_row=[]
        run_accession_new_row.append(df['run_accession'].iloc[i].tolist())
        run_accession_new_row=list(chain.from_iterable(run_accession_new_row))
        run_accession_new.append(run_accession_new_row)
        #this also needs the read_files

        read_file_new_row=[]
        read_file_new_row.append(df['read_file'].iloc[i].tolist())
        read_file_new_row=list(chain.from_iterable(read_file_new_row))
        read_file_new.append(read_file_new_row)


    replicates_list=list(range(0, replicate_nr))

    replicates_columns=[]
    for i in replicates_list:
        replicates_columns.append(str(i))
    print(replicates_columns)
    naming=["sample_accession_repl", "run_accession_repl","read_file_repl"] #3 three times
    out=list(itertools.product(naming, replicates_columns))

    detouple=()
    for i in out:
        detouple+=i
    detouple=iter(detouple)
    print(detouple)
    replicates_columns=[i+next(detouple, '') for i in detouple]


    sample_accession_new=pd.DataFrame(sample_accession_new)

    run_accession_new=pd.DataFrame(run_accession_new)

    read_file_new=pd.DataFrame(read_file_new)

    print(read_file_new)


        #merge all together unique/dropna and the between newly created ones
    out_df=pd.concat([sample_accession_new, run_accession_new, read_file_new], axis=1)


    print(out_df)
    print(replicates_columns)
    out_df.columns=replicates_columns
    out_df["sample_accession"]=out_df.filter(regex=("sample_accession*")).values.tolist()
    out_df["run_accession"]=out_df.filter(regex=("run_accession*")).values.tolist()
    out_df["read_file"]=out_df.filter(regex=("read_file*")).values.tolist()
    #out_df.filter(regex=("d.*"))
    out_df.drop(replicates_columns, axis=1, inplace=True)
    #merge all together unique/dropna and the two newly created ones


    df=df.drop(keep_list, axis=1)

    df=df.drop_duplicates(keep='first', inplace=False)
    df=df.reset_index(drop=True)

    out_df=pd.concat([df, out_df], axis=1)

    return(out_df)


def drop_redundant_v2(df, redundant_groups):


    redundant_groups = [item for sublist in redundant_groups for item in sublist] #flattens
    redundant_groups.sort()



    if not redundant_groups:
    #    print("nothing in redundant list none dropped")
        return(df)
    elif len(df.columns)==1:
        return(df)
    else:

        #print(list_)
        #red_gr = [item for sublist in list_ for item in sublist]
        #red_gr=list(chain.from_iterable(list))

        df = df.reindex(sorted(df.columns), axis=1)
        #print("old_df")
    #    print(df)
        redundant_groups=df.columns.tolist()
        #red_gr.sort()
        #print(red_gr)

        to_be_merged_df=df #this I guess is unncessary as it does not affect the size while grouping
        #print(list_)
        #print(redundant_groups)
        red_pairs=itertools.combinations(redundant_groups, r=2)

        #print(red_pairs)

        red_pairs_list=[]
        for i in red_pairs:
            temp=tuple(i)
            red_pairs_list.append(temp)
        #print(red_pairs_list)

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

            elif len(to_be_merged_df.groupby(f'{i[0]}')) == len(to_be_merged_df.groupby([f'{i[0]}', f'{i[1]}'])):

                 to_be_merged_df.drop([f'{i[0]}'], axis=1, inplace=True)


                 return(df)

            else:

                return(df)

def detect_multiple_studies(df):
    for y, x in df.groupby(df.isnull().dot(df.columns)):


        substudies = {x : y.dropna(1) for x, y in df.groupby(df.isnull().dot(df.columns))}


    print("Substudies found:")
    substudies_names=list(reversed(list(substudies.keys())))
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
        #print(level_size)
        replicates_nr=replicates_def(sample_attributes, level_size)
    else:

        factors_levels=factor_levels_dict(sample_attributes)
        level_size=factors_levels_size(factors_levels)

        replicates_nr=replicates_def(sample_attributes, level_size)
    if replicates_nr[0]: #true there are replicates
        sample_attributes=replicates_merger(sample_attributes, list(factors_levels.keys()), replicates_nr[1])
    else:
        sample_attributes=sample_attributes


    print(f"{study} ")
    read_type_det
    global paired
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

    print(f"Number of Replicates: {replicates_nr[1]}.")
    print(f"Experimental factors are {list(factors_levels.keys())}")
    print(f"with levels {list(factors_levels.values())}")
    print("\nSample table:")

    print(sample_attributes[sample_attributes.columns[~sample_attributes.columns.isin(sample_info_col_list)]])
    print(f"\nPipeline will investigate Differential Gene/Transcript Expression (DGE/DTE) and \nnovel and annotated Alternative Splicing (AS) events\nusing pairwise comparisons between all conditions.\n")
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

            confirm_subset_by_level=input(f"Study will be subsetted by level {subset_by_level} enter y to continue, n to abort: ")
            if confirm_subset_by_level=="y":

                subsetted_sample_attributes=subset_sample_table(subset_by_level, sample_attributes)

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

        sample_attributes_confirmed=sample_attributes
        contrast_indices_confirmed=contrasts_indices

        sample_attributes_confirmed.reset_index(drop=True, inplace=True)


        comp_A=[]
        comp_B=[]


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
        comp_A_names=[i.replace(" ", "") for i in comp_A_names]
        comp_A_names=[i.replace("/", "_") for i in comp_A_names]
        comp_A_names=[i.replace("(", "_") for i in comp_A_names]
        comp_A_names=[i.replace(")", "_") for i in comp_A_names]
        comp_A_names = [i.encode('ascii', 'ignore').decode('ascii') for i in comp_A_names] #unicode to byte  conversion is a problem for pandas
        print("First pair of comparison: ", comp_A_names)
        print("\n")
        run_list_A=[run_list[i] for i in A]
        B=[i[1] for i in contrast_indices_confirmed]
        comp_B_names=[comp[i] for i in B]
        comp_B_names=[i.replace(" ", "") for i in comp_B_names]
        comp_B_names=[i.replace("/", "_") for i in comp_B_names]
        comp_B_names=[i.replace("(", "_") for i in comp_B_names]
        comp_B_names=[i.replace(")", "_") for i in comp_B_names]
        comp_B_names = [i.encode('ascii', 'ignore').decode('ascii') for i in comp_B_names]
        print("Second pair of comparison:", comp_B_names)
        run_list_B=[run_list[i] for i in B]

        sample_acc_list_A=[sample_acc_list[i] for i in A]
        sample_acc_list_B=[sample_acc_list[i] for i in B]


        run_accession_list_A=[run_accession_list[i] for i in A]
        run_accession_list_B=[run_accession_list[i] for i in B]


        ################ names ############
        contrast_names = [i + "_vs_"+ j for i, j in zip(comp_A_names, comp_B_names)]
        contrast_names=[i.replace(" ", "") for i in contrast_names]
        contrast_names=[i.replace("/", "_") for i in contrast_names]
        #print(contrast_names)

        ############# bam dict ############
        print("bam_dict")
        print(replicates_nr[1])
        print(range(0,replicates_nr[1]))
        replicate_suffix=list(range(0,replicates_nr[1]))*len(comp_A_names)
        print("suff", replicate_suffix)
        comp_A_names_r=[item for item in comp_A_names for i in range(0,replicates_nr[1])]
        print("reps", comp_A_names_r)
        comp_A_names_replicates=[i+ "_R" + str(j) for i,j in zip(comp_A_names_r,replicate_suffix)]

        comp_A_names_replicates = [comp_A_names_replicates[i:i+replicates_nr[1]] for i in range(0, len(comp_A_names_replicates), replicates_nr[1])]
        print(comp_A_names_replicates)

        run_accession_list_A
        cond_A_acc=zip(comp_A_names_replicates, run_accession_list_A)
        print("zipped", cond_A_acc)
        comp_B_names_r=[item for item in comp_B_names for i in range(0,replicates_nr[1])]
        comp_B_names_replicates=[i+ "_R" + str(j) for i,j in zip(comp_B_names_r,replicate_suffix)]
        comp_B_names_replicates = [comp_B_names_replicates[i:i+replicates_nr[1]] for i in range(0, len(comp_B_names_replicates), replicates_nr[1])]

        run_accession_list_B
        cond_B_acc=zip(comp_B_names_replicates, run_accession_list_B)

        key_1 = 'A'
        key_2 = 'B'
        comparison_dict_v2={k : {key_1 : v1, key_2 : v2} for k,v1,v2 in zip(contrast_names, cond_A_acc, cond_B_acc)}
        print(comparison_dict_v2)


        A_list=list(zip(comp_A_names, sample_acc_list_A, run_list_A))

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

        #################

        key_1 = 'A'
        key_2 = 'B'
        comparison_dict={k : {key_1 : v1, key_2 : v2} for k,v1,v2 in zip(contrast_names, run_accession_list_A, run_accession_list_B)}




        print("Number of comparisons:", len(contrasts_indices))
        print("\n")
        print(replicates_nr[0])
        start_pipeline=input("To issue pipleline using pairwise comparison between each to HPC enter y, to abort, enter n: ")
        if start_pipeline=="y":
            study_total=open(output_dir +'study_total.pkl', 'wb')
            replicates=replicates_nr[0]
            pickle.dump(study, study_total)
            pickle.dump(tax_id, study_total)
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
            print("Input saved.")
            os.system("sbatch snake.sh")
            return(sample_attributes, contrasts_indices)



            not_started=False
        if start_pipeline=="n":
            exit()

if args.ena:
  study=args.ena #str(sys.argv[1])
  study=''.join(study)
  cwd = os.getcwd()

  output_dir=cwd+"/"
  if not os.path.exists(output_dir):
      os.makedirs(output_dir)

  print(f"Retrieving information for study {study} from ENA")


  pd.set_option('display.max_columns', 500)
  pd.set_option('display.width', 1000)
  study_table_request = requests.get(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={study}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true")
  print(study_table_request)
  study_table=study_table_request.content

  study_table=str(study_table, 'utf-8')
  print(study_table)
  study_table=StringIO(study_table)
  print(study_table)
  study_table=pd.read_csv(study_table, sep='\t', lineterminator='\n')
  print(study_table)
  tax_id=study_table['tax_id'].tolist()
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

          first=i[0]
          first=os.path.basename(first)
          read_type.append(first)






  run_accession=study_table['run_accession'].tolist()

  xml_links=[]
  for i in sample_accession:
      xml_links.append(f"https://www.ebi.ac.uk/ena/browser/api/xml/{i}")



  xml_list= [urlopen(link).read() for link in xml_links]


  xtree=[et.ElementTree(et.fromstring(xml)) for xml in xml_list]


  xroot=[xtree_i.getroot() for xtree_i in xtree]


  sample_attributes=[]
  for i in xroot:
      xheader_inner=[]
      xvalues_inner=[]
      for j in i.findall(".//SAMPLE_ATTRIBUTE/TAG"):

          xheader_inner.append(j.text)
      for j in i.findall(".//SAMPLE_ATTRIBUTE/VALUE"):

          xvalues_inner.append(j.text)
      attributes_inner=dict(zip(xheader_inner, xvalues_inner))
      sample_attributes.append(attributes_inner)

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
  print("b4")
  print(sample_attributes)
  sample_attributes=pd.DataFrame(sample_attributes)


  sample_attributes = sample_attributes[sample_attributes.columns.drop(list(sample_attributes.filter(regex='ENA'))) ]
  print("here")
  print(sample_attributes)

  sample_info_col_list=['sample_accession', 'run_accession', 'read_file']

  columns_useable=sample_attributes.columns[~sample_attributes.columns.isin(sample_info_col_list)]

  print(columns_useable)

  for col in columns_useable:

      if len(sample_attributes[col].unique()) == 1:
          sample_attributes.drop(col,inplace=True,axis=1)




  sample_attributes_confirmed, contrast_indices_confirmed=contrast_prompter()

if args.local:
   print("local mode")
   print(''.join(args.local))

   study_table=pd.read_csv(''.join(args.local), sep=',', lineterminator='\n')
   print(study_table)
   study=study_table['study_accession'].tolist()[0]
   tax_id=study_table['tax_id'].tolist()
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
           first=i[0]
           first=os.path.basename(first)
           read_type.append(first)






   run_accession=study_table['run_accession'].tolist()
   sample_attributes=study_table
   sample_info_col_list=['sample_accession', 'run_accession', 'read_file']

   columns_useable=sample_attributes.columns[~sample_attributes.columns.isin(sample_info_col_list)]



   for col in columns_useable:

      if len(sample_attributes[col].unique()) == 1:
         sample_attributes.drop(col,inplace=True,axis=1)


   cwd = os.getcwd()

   output_dir=cwd+"/"
   if not os.path.exists(output_dir):
      os.makedirs(output_dir)


   sample_attributes_confirmed, contrast_indices_confirmed=contrast_prompter()
