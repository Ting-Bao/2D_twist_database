# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com

import os
import json
from pymatgen.io.vasp import Poscar


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

def check_ehull(jsonpath,ehull=0.1):
    '''read xx.alldata.json, check energy above hull, should be <0.1eV/atom
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

def check_hexagonal(poscarpath):
    '''check if the poscar is hexagional structure
    https://zhuanlan.zhihu.com/p/292407444
    to distinguish the alpha 60/120 degree, angle alpha is also returned
    return (bool, angle alpha of poscar)
    '''
    struc = Poscar.from_file(poscarpath)
    temp = struc.structure.lattice
    return temp.is_hexagonal,temp.angles

def check_square(poscarpath):
    struc = Poscar.from_file(poscarpath)
    temp = struc.structure.lattice
    return temp.is_orthogonal and abs(temp.a-temp.b)<0.001

def runsh_w001(name='default'):
    with open('template/VASP/run_vasp_w001.sh') as f:
        temp=f.readlines()
    for i in range(len(temp)):
        temp[i]=temp[i].replace('NAMENAME', name)
    return temp
