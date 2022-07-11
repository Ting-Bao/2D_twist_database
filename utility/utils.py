# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com

import os
import json
import re
import pandas as pd
from pymatgen.io.vasp import Poscar, Outcar
from monty.io import reverse_readfile

def check_path(path,creat_if_not=True):
    if not os.path.exists(path):
        if creat_if_not:
            os.makedirs(path)
        return False
    return True

def read_json(jsonpath):
    '''return the content of a jsonfile
    '''
    with open(jsonpath,'r',encoding='utf-8') as f:
        temp=json.load(f)
    return temp
    
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

def opmxbatch_w001(name='default'):
    with open('template/VASP/run_vasp_w001.sh') as f:
        temp=f.readlines()
    for i in range(len(temp)):
        temp[i]=temp[i].replace('NAMENAME', name)
    return temp

def find_atomic_number(poscarpath):
    '''return a list containing atomic numbers
    '''
    struc = Poscar.from_file(poscarpath)
    temp = struc.structure.atomic_numbers

    return list(temp)

def grep_TOTEN(outcarpath):
    ''' get TOTEN from outcar
    return the value, if not found, return None
    faster than using pymatgen's outcar.final_energy which requires time-consuming initialization
    the code here refer to pymatgen's source code
    '''
    temp=[]
    e_fr_energy_pattern=re.compile(r"free  energy   TOTEN\s+=\s+([\d\-\.]+)")
    for line in reverse_readfile(outcarpath):
        clean = line.strip() 
        m = e_fr_energy_pattern.search(clean)
        if m:
            temp.append(m.group(1))
        if len(temp)>0:
            return float(temp[0])
    return None
    # following code is too time consuming due to class Outcar's initialization
            #outcar=Outcar(filename+'/OUTCAR')
            #print(outcar.final_energy)

def read_excel(excelpath):
    '''read excel (.xlsx) and return the pd.DataFrame
    '''
    temp=pd.read_excel(excelpath,index_col=0)
    return temp

def add_to_json(jsonpath,key,content):
    '''the key should have a list type value in the json file
    this func append the content to the correspond list
    '''
    with open(jsonpath,'r') as f:
        temp=json.load(f)
        temp[key].append(content)
    with open(jsonpath,'w',encoding='utf-8') as f:
        json.dump(temp,f,indent=2)
    return temp[key]