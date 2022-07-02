# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# 
# here we decide automatically if a material should consider SOC during SSC calc
# rules are as follows:
# 1. if all elements contained has atomic number <=18, then donot apply soc
# 2. if both homo and lumo diff < 0.02, then donot apply soc 
# 
# !!! human experts can further judge the results !!!

from utility.utils import *
import os
import numpy as np
import json
import shutil
import pandas as pd

used_data='data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
socfilepath='data/ifsoc_data/'

tofile='data/ifsoc_data/AAresult/'
finalsocfile='data/finalchoice/ifsoc.json'

plotfile='utility/plot_openmx_band.py'
comparefile='utility/compare_band.py'

def rule1(poscarpath):
    # return if all elements contained has atomic number < 18
    temp=find_atomic_number(poscarpath)
    if max(temp)<19:
        return False
    return True

def rule2(comparejson):
    with open(comparejson, 'r') as f:
        temp=json.load(f)
    if temp['lumo_diff']<0.02 and temp['homo_diff']<0.02:
        return False # false means not need for soc
    return True

def main():
    check_path(tofile)
    namelist = read_namelist(used_data+'namelist.json')
    
    # a dict to store if a material should apply soc
    socresult={}

    #gen openmx.Band.json and copy to new file to compare them 
    for i in namelist:
        check_path(tofile+i+'/')
        if not os.path.exists(socfilepath+i+'_soc/openmx.Band.json'):
            os.system('cd {} && python {}'.format(socfilepath+i+'_soc','../../../'+plotfile))
        if not os.path.exists(socfilepath+i+'_nosoc/openmx.Band.json'):
            os.system('cd {} && python {}'.format(socfilepath+i+'_nosoc','../../../'+plotfile))
        if not os.path.exists(tofile+i+'/nosoc-soc.compare.json'):
            shutil.copy(socfilepath+i+'_soc/openmx.Band.json',tofile+i+'/soc.Band.json')
            shutil.copy(socfilepath+i+'_nosoc/openmx.Band.json',tofile+i+'/nosoc.Band.json')
            os.system('cd {}&& python {} -u 6.0 -d -6.0'.format(tofile+i,'../../../../'+comparefile))  

    for i in namelist:
        if rule1(socfilepath+i+'_soc/POSCAR')==False:
            socresult[i]=False
        else:
            socresult[i]=rule2(tofile+i+'/nosoc-soc.compare.json')

    pd.DataFrame.from_dict(socresult,orient='index').to_excel(tofile+'data.xlsx')
    pd.DataFrame.from_dict(socresult,orient='index').to_excel('data/finalchoice/'+'ifsoc.xlsx')
    with open(finalsocfile,'w',encoding='utf-8') as f:
        json.dump(socresult,f)

if __name__=='__main__':
    main()