# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# 
# here we get the most favorable stacking of materials
# require calculated stacking file using  VASP

# we also pick the priority due to ICSD id
# put the final choice into 2 sets: with/without ICSD id


from utility.utils import *
import os
import numpy as np
import pandas as pd
import shutil
import json


used_data='data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
filepath='data/calculated_ABstacking_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
alldatapath='data/c2dbdata/jsondata/'
socfile='data/finalchoice/ifsoc.json'

topath='data/finalchoice/all/'
# where to put the final choice of the materials

typelist=['AA','AB','AAp1','ABp1','AAp2','ABp2'] # put from most favorable to most unfavorable 
compare_tol=0.0 # unit:eV
# see try_stacking.py to find the definition 

result_excel='data/finalchoice/stacking_choice.xlsx'

def main():
    check_path(topath)
    namelist=read_namelist(used_data+'namelist.json')
    #print(namelist)
    Etable=pd.DataFrame(columns=typelist+['favorate_stack'])

    # check the TOTEN, collect into a table
    # check the most favorable stacking, a tolerence is applied
    for i in namelist:
        temp=()
        lowest=0
        for j in range(len(typelist)):
            filename=filepath+i+'_'+typelist[j]
            TOTEN=grep_TOTEN(filename+'/OUTCAR')
            if TOTEN==None:
                #here means the calculation dosen't normally finished
                print(i,typelist[j],TOTEN)
            temp=temp+(TOTEN,)
            if TOTEN!=None:
                if temp[lowest]==None or temp[lowest]-temp[j]>compare_tol:
                    lowest=j
        item=pd.DataFrame([temp+(typelist[lowest],)],columns=typelist+['favorate_stack'],index=[i])
        Etable=pd.concat([Etable,item],ignore_index=False)
        # to be done: get soc and put into the table
    Etable.to_excel(result_excel)
    print(Etable)


    ####################
    # pick the priority
    ####################
    #read the soc info,use human 
    socinfo=read_json(socfile)
    priority_list=[]
    normal_list=[]
    path_p=topath+'../part1_prioirty/'
    path_n=topath+'../part2_normal/'
    check_path(path_p)
    check_path(path_n)
    
    for i in namelist:
        choice = Etable.at[i,'favorate_stack']
        choicefile = filepath+i+'_'+choice
        check_path(topath+i)
        if not os.path.exists(topath+i+'/relaxed_POSCAR'):
            shutil.copy(choicefile+'/CONTCAR', topath+i+'/relaxed_POSCAR')
        if not os.path.exists(topath+i+'/prepare.json'):
            prepare_info={'from_stack':choice}
            prepare_info['soc']=socinfo[i]
            with open(topath+i+'/prepare.json', 'w',encoding='utf-8') as f:
                json.dump(prepare_info,f)
        # put into priority/normal list according to ICSD id
        if if_ICSD(alldatapath+i+'.all_data.json'):
            priority_list.append(i)
            shutil.copytree(topath+i,path_p+i)
        else:
            normal_list.append(i)
            shutil.copytree(topath+i,path_n+i)


if __name__=='__main__':
    main()