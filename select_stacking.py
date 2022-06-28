# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# 
# here we get the most favorable stacking of materials
# require calculated stacking file using  VASP

from utility.utils import *
import os
import numpy as np

used_data='data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
filepath='data/ABstacking_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
alldatapath='data/c2dbdata/jsondata/'

topath='data/finalchoice/all/'
# where to put the final choice of the materials

type=['AA','AAp1','AAp2','AB','ABp1','ABp2']
# see try_stacking.py to find the definition 

if __name__=='__main__':
    check_path(topath)
    namelist=read_namelist(used_data+'namelist.json')
    print(namelist)


