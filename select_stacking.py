# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# 
# here we get the most favorable stacking of materials
# require calculated stacking file using  VASP

from utility.utils import *
import os
import numpy as np
import pandas as pd
import shutil

used_data='data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
filepath='data/calculated_ABstacking_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
alldatapath='data/c2dbdata/jsondata/'

topath='data/finalchoice/all/'
# where to put the final choice of the materials

typelist=['AA','AB','AAp1','ABp1','AAp2','ABp2'] # put from most favorable to most unfavorable 
compare_tol=0.0 # unit:eV
# see try_stacking.py to find the definition 


if __name__=='__main__':
    check_path(topath)
    namelist=read_namelist(used_data+'namelist.json')
    #print(namelist)
    Etable=pd.DataFrame(columns=typelist+['favorate_stack'])

    # check the TOTEN, collect into a table
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
    Etable.to_excel('data/finalchoice/stacking_choice.xlsx')
    print(Etable)
    

    # check the most favorable stacking, a tolerence is applied
    #for i in namelist:
    #    print('')
            

            

