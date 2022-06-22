# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
#
# use the primtive cell of selected data, use openmx to get the band structure with/without soc
# compare the band to decide whether to consider SOC in the dataset generattion


from utility.utils import *
import shutil
import os

loc=os.getcwd()
used_data=loc+'/data/selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
topath=loc+'/data/ifsoc_data/'
gencode=loc+'/utility/gen_openmx_input_dftu_v5x.py'


if __name__=='__main__':
    namelist=read_namelist(used_data+'namelist.json')
    # print(namelist,len(namelist))
    for i in range(4):
        name=namelist[i]
        print (name)
        if not os.path.exists(topath+name+'_soc'):
            os.makedirs(topath+name+'_soc')
        if not os.path.exists(topath+name+'_nosoc'):
            os.makedirs(topath+name+'_nosoc')
        shutil.copy(used_data+'poscar/POSCAR_{}'.format(name), topath+name+'_soc/POSCAR')
        shutil.copy(used_data+'poscar/POSCAR_{}'.format(name), topath+name+'_nosoc/POSCAR')
        os.system('cd {}_soc && python {} --soc True'.format(topath+name,gencode))
        os.system('cd {}_nosoc && python {} --soc False'.format(topath+name,gencode))
        shutil.copy('template/OPMX/run_opmx_w001.sh', topath+name+'_soc')
        shutil.copy('template/OPMX/run_opmx_w001.sh', topath+name+'_nosoc')