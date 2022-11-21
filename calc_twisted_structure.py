# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.11
#
# this is to generate the 1x2 and part of 2x3 structures and conduct openmx calculation

namelist='./data/finalchoice/all/namelist.json'
frompath='./data/finalchoice/all/'
topath='./data/twisted_band_opmx/'


from utility.utils import *
from utility.gen_opmx_input_func import *
from utility.gen_twist.gen_twist import *
import os
import shutil
import json

namelist=read_namelist(namelist)
#print(namelist,len(namelist))

for name in namelist:
    if not os.path.exists(topath+'{}_2-1'.format(name)):
        os.makedirs(topath+'{}_2-1'.format(name))
        twisted_angles, atom_num = gen_twist(m=2, n=1, fromfile='data/finalchoice/all/{}/relaxed_POSCAR'.format(name),
                tofile=topath+'{}_2-1/POSCAR'.format(name))
        shutil.copy('data/finalchoice/all/{}/prepare.json'.format(name),topath+'{}_2-1/prepare.json'.format(name))
        shutil.copy('template/OPMX/run_opmx_w001.sh'.format(name),topath+'{}_2-1/'.format(name))

    if not os.path.exists(topath+'{}_3-2'.format(name)):
        os.makedirs(topath+'{}_3-2'.format(name))
        twisted_angles, atom_num = gen_twist(m=3, n=2, fromfile='data/finalchoice/all/{}/relaxed_POSCAR'.format(name),
                tofile=topath+'{}_3-2/POSCAR'.format(name))
        shutil.copy('data/finalchoice/all/{}/prepare.json'.format(name),topath+'{}_3-2/prepare.json'.format(name))
        shutil.copy('template/OPMX/run_opmx_w001.sh'.format(name),topath+'{}_3-2/'.format(name))
        if atom_num>120:
            os.system('rm -r '+ topath +'{}_3-2/'.format(name))

filelist=[]
for i in os.listdir(topath):
    if os.path.isdir(topath+i):
            filelist.append(topath+i)

#print(filelist)
runall=[]
absolute_path='/home/xurz/temp_baot/twisted_band_opmx/'
for i in filelist:
    with open(i+'/prepare.json','r',encoding='utf-8') as f:
        temp=json.load(f)
    soctag=temp['soc']
    gen_openmx_input(frompath=i+'/POSCAR',topath=i+'/',soc=soctag)
    runall.append('cd {} && qsub run_opmx_w001.sh\n'.format(absolute_path+i))

with open(topath+'runall.sh', 'w', encoding = 'utf-8') as f:
    f.writelines(runall)


