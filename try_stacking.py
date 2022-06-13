# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
#
# for hexagonal lattice try AA and AB stacking and find the lower one
# VASP calculation is required here, thus the input file is generated locally and should be uploaded to the server for calculation
# vaspkit is required to get POTCAR and KPOINTS, make sure vaspkit have proper access to PBE pseudo-potential

used_data='selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
vaspkitpath='/home/tingbao/Software/vaspkit.1.3.3/bin/vaspkit' 
topath='ABstacking_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'

import os
import shutil
import json
from make_bilayer_v2 import *

def read_namelist(namelistpath):
    with open(namelistpath,'r') as f:
        temp=json.load(f)
    return temp['namelist']

if __name__=='__main__':
    os.makedirs(topath+'AB/')
    os.makedirs(topath+'AA/')

    namelist=read_namelist(used_data+'namelist.json')
    for i in range(len(namelist)):
        name=namelist[i]
        poscarpath=used_data+'poscar/POSCAR_{}'.format(name)
