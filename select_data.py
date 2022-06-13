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


poscarpath='c2dbdata/poscardata/'
alldatapath='c2dbdata/jsondata/'
max_ele_type=2
max_atoms=8
check_ICSD=False
ehull=0.1

import os
import json
from pymatgen.io.vasp import Poscar
from tqdm import tqdm
import shutil

def read_namelist(namelistpath):
    with open(namelistpath,'r') as f:
        temp=json.load(f)
    return temp['namelist']

def if_ICSD(jsonpath):
    '''judge whether a material has a ICSD id, indicating its existance in real wold.'''
    with open(jsonpath,'r') as f:
        temp=json.load(f)
    if "icsd_id" in temp['info.json'].keys():
        return True
    return False

def check_ehull(jsonpath,ehull=ehull):
    '''check energy above hull, should be <0.1eV/atom
    '''
    with open(jsonpath,'r') as f:
        temp=json.load(f)
    #incase no ehull tag
    try:
        if temp['results-asr.convex_hull.json']['kwargs']['data']['ehull']<ehull:
            return True
        return False
    except:
        return False

def count_atom(poscarpath):
    '''return number of atoms in the poscar'''
    struc = Poscar.from_file(poscarpath)
    return len(struc.as_dict()['structure']['sites'])

def count_ele(poscarpath):
    '''return number of types of element in the poscar.'''
    struc = Poscar.from_file(poscarpath)
    ele=[]
    for i in struc.as_dict()['structure']['sites']:
        temp=i['species'][0]['element']
        if temp not in ele:
            ele.append(temp)
    return len(ele)

def check_hexagonal(poscar):
    '''check if the poscar is hexagional structure
    https://zhuanlan.zhihu.com/p/292407444
    to distinguish the alpha 60/120 degree, angle alpha is also returned
    return (bool, angle alpha of poscar)
    '''
    struc = Poscar.from_file(poscarpath)
    temp = struc.structure.lattice
    return temp.is_hexagonal(),temp.angles


if __name__=='__main__':
    #where to put the selected data
    configname='ICSD_{}+maxele_{}+maxatom_{}+ehull_{:.1f}'.format(check_ICSD,max_ele_type,max_atoms,ehull)
    jsdst='selected_data/{}/jsonfile/'.format(configname)
    psdst='selected_data/{}/poscar/'.format(configname)
    os.makedirs(jsdst)
    os.makedirs(psdst)

    namelist=read_namelist('c2dbdata/namelist.json')
    print ('we have {} candidates.'.format(len(namelist)))

    count=0
    newnamelist=[]
    for i in range(len(namelist)):
        name=namelist[i]
        jsonpath='c2dbdata/jsondata/{}.all_data.json'.format(name)
        poscarpath='c2dbdata/poscardata/POSCAR_{}'.format(name)
        if check_ICSD and not if_ICSD(jsonpath):
            print (name,'\t no ICSD ID, passed')
            continue 
        if not check_ehull(jsonpath):  
            print (name,'\t energy above hull > 0.1, passed')
            continue  
        if count_ele(poscarpath) > max_ele_type:
            print (name,'\t too many element types, passed')
            continue
        if count_atom(poscarpath) > max_atoms:
            print (name,'\t too many atoms, passed')
            continue
        hexagonal=check_hexagonal(poscarpath)
        if hexagonal[0]:
            print (hexagonal)
        else:
            print (name,' not hexagonal, passed')
            continue
        print (name, '\t satisfies all conditions, OK!')
        newnamelist.append(name)
        count+=1
        shutil.copy(jsonpath.format(name),jsdst)
        shutil.copy(poscarpath, psdst)
    with open('selected_data/{}/namelist.json'.format(configname),'w',encoding='utf-8') as f:
        json.dump({'namelist':newnamelist},f,indent=2)
    print('Finally got {} materials, copied to selected data'.format(count))