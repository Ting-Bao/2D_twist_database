# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
# 
# use this code to find if a material requires SOC adjustment in electronic structrue calcualtion
# use the primtive cell of selected data, use openmx to get the band structure with/without soc
# compare the band to decide whether to consider SOC in the dataset generattion


from utility.utils import *
from utility.gen_opmx_input_func import *
import shutil
import os

loc=os.getcwd()
used_data=loc+'/data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
topath=loc+'/data/ifsoc_data/'
gencode=loc+'/utility/gen_openmx_input_dftu_v6_local.py'


def main():
    namelist=read_namelist(used_data+'namelist.json')
    # print(namelist,len(namelist))
    for i in range(len(namelist)):
        name=namelist[i]
        print (name)
        if not os.path.exists(topath+name+'_soc'):
            os.makedirs(topath+name+'_soc')
        if not os.path.exists(topath+name+'_nosoc'):
            os.makedirs(topath+name+'_nosoc')
        shutil.copy(used_data+'poscar/POSCAR_{}'.format(name), topath+name+'_soc/POSCAR')
        shutil.copy(used_data+'poscar/POSCAR_{}'.format(name), topath+name+'_nosoc/POSCAR')
        gen_openmx_input(frompath=topath+name+'_nosoc/POSCAR',topath=topath+name+'_nosoc/',soc=False)
        gen_openmx_input(frompath=topath+name+'_soc/POSCAR',topath=topath+name+'_soc/',soc=True)
        shutil.copy('template/OPMX/run_opmx_w001.sh', topath+name+'_soc')
        shutil.copy('template/OPMX/run_opmx_w001.sh', topath+name+'_nosoc')
    
    with open(topath+'runall.sh','w') as f:
        temp=[]
        for i in namelist:
            temp.append('cd {}_nosoc && qsub run_opmx_w001.sh &&cd ..\n'.format(i))
            temp.append('cd {}_soc && qsub run_opmx_w001.sh &&cd ..\n'.format(i))
        f.writelines(temp)

if __name__=='__main__':
    main()