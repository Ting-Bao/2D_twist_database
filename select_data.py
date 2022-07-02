# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
#
# require the existed c2db database with rough selection of stable and nonmagnetic already done
# which have xx.alldata.json and POSCAR_xx, where xx are corresponded
#
# here, we select further by:
# 1. should have ICSD ID
# 2. should have elements type less than N (N=2)
# should be hexagonal lattice, which means:
#       3. should have lattice vector |a|==|b|
#       4. currently we only take gamma=120/60 degree


poscarpath='data/c2dbdata/poscardata/'
alldatapath='data/c2dbdata/jsondata/'
max_ele_type=2
max_atoms=8
check_ICSD=False
ehull=0.1

import os
import json
from pymatgen.io.vasp import Poscar
from tqdm import tqdm
from utility.utils import *
import shutil


def main():
    #where to put the selected data
    configname='square_ICSD_{}+maxele_{}+maxatom_{}+ehull_{:.1f}'.format(check_ICSD,max_ele_type,max_atoms,ehull)
    jsdst='data/selected_data/{}/jsonfile/'.format(configname)
    psdst='data/selected_data/{}/poscar/'.format(configname)
    os.makedirs(jsdst)
    os.makedirs(psdst)

    namelist=read_namelist('data/c2dbdata/namelist.json')
    print ('we have {} candidates.'.format(len(namelist)))

    count=0
    newnamelist=[]
    for i in range(len(namelist)):
        name=namelist[i]
        jsonpath='data/c2dbdata/jsondata/{}.all_data.json'.format(name)
        poscarpath='data/c2dbdata/poscardata/POSCAR_{}'.format(name)
        if check_ICSD and not if_ICSD(jsonpath):
            print (name,'\t no ICSD ID, out')
            continue 
        if not check_ehull(jsonpath,ehull=ehull):  
            print (name,'\t energy above hull > 0.1, out')
            continue  
        if count_ele(poscarpath) > max_ele_type:
            print (name,'\t too many element types, out')
            continue
        if count_atom(poscarpath) > max_atoms:
            print (name,'\t too many atoms, out')
            continue
        '''
        hexagonal=check_hexagonal(poscarpath)
        if hexagonal[0]:
            print (hexagonal)
        else:
            print (name,' not hexagonal, out')
            continue
        '''
        if not check_square(poscarpath):
            continue
        print (name, '\t satisfies all conditions, OK!')
        newnamelist.append(name)
        count+=1
        shutil.copy(jsonpath.format(name),jsdst)
        shutil.copy(poscarpath, psdst)
    with open('data/selected_data/{}/namelist.json'.format(configname),'w',encoding='utf-8') as f:
        json.dump({'namelist':newnamelist},f,indent=2)
    print('Finally got {} materials, copied to selected data'.format(count))

if __name__=='__main__':
    main()