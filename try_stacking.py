# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
#
# for hexagonal lattice try following AA,AAp,AB... stacking and find the one with lower energy:
# here angle gamma = 120 degree
# AA: upper == lowver
# AAp1: upper == lower + (0.33a,0.67b) translation
# AAp2: upper == lower + (0.67a,0.33b) translation
# AB: upper == lower rotate 180 degree
# ABp1: upper == AB lower + (0.33a,0.67b) translation
# ABp2: upper == AB lower + (0.67a,0.33b) translation
# 
# VASP calculation is required here, thus the input file is generated locally and should be uploaded to the server for calculation
# vaspkit is required to get POTCAR and KPOINTS, make sure vaspkit have proper access to PBE pseudo-potential

used_data='selected_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'
vaspkitpath='/home/tingbao/Software/vaspkit.1.3.3/bin/vaspkit' 
savepath='ABstacking_data/ICSD_False+maxele_3+maxatom_8+ehull_0.1/'

import os
import shutil
import json
from make_bilayer_v2 import *
from utils import *

def prepare_for_VASP(calcpath,name='default'):
    '''POSCAR file should be in the calcpath
    '''
    run_sh=runsh_w001(name=name)
    # copy incar to calc folder
    shutil.copy('template/VASP/INCAR.RELAX', toposcar_aa+'INCAR')
    shutil.copy('template/VASP/INCAR.RELAX', toposcar_ab+'INCAR')
    # use vaspkit to generate potcar and kpoints
    os.system('cd {} &&  {} -task 103'.format(calcpath,vaspkitpath)) # POTCAR
    os.system('cd {} &&  (echo 102; echo 2;echo 0.02)|{} '.format(calcpath,vaspkitpath)) #KPOINTS
    with open(calcpath+'run.sh','w',encoding='utf-8') as f:
        f.writelines(run_sh)
    
    print("!!!!!!!!!!{}!!!!!!!!!!!!".format(calcpath))
    # in vaspkit 1.3, this command dosen't work: vaspkit -task 102 -kps Gamma -kpr 0.02
    #print (calcpath, ' prepared')

if __name__=='__main__':
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    
    namelist=read_namelist(used_data+'namelist.json')

    for i in range(len(namelist)):
        name=namelist[i]
        fromposcar=used_data+'poscar/POSCAR_{}'.format(name)
        toposcar_aa=savepath+name+'_AA/'
        toposcar_ab=savepath+name+'_AB/'
        toposcar_aap1=savepath+'AAAtemp/'+name+'_AAp1/'
        toposcar_aap2=savepath+'AAAtemp/'+name+'_AAp2/'
        toposcar_abp1=savepath+'AAAtemp/'+name+'_ABp1/'
        toposcar_abp2=savepath+'AAAtemp/'+name+'_ABp2/'

        toposcarlist=[toposcar_aa,toposcar_ab,toposcar_aap1,toposcar_aap2,toposcar_abp1,toposcar_abp2]

        for i in toposcarlist:
            if not os.path.exists(i):
                os.makedirs(i)
        
        # make the stacking poscars
        # AA AB AAp1 AAp2 ABp1 ABp2
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_aa+'POSCAR',rotate_angle=0)
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_ab+'POSCAR',rotate_angle=180)
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_aap1+'POSCAR',ashift=0.3333,bshift=0.6667,rotate_angle=0)
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_aap2+'POSCAR',ashift=0.6667,bshift=0.3333,rotate_angle=0)
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_abp1+'POSCAR',ashift=0.3333,bshift=0.6667,rotate_angle=180)
        make_bilayer(fromposcar=fromposcar, toposcar=toposcar_abp2+'POSCAR',ashift=0.6667,bshift=0.3333,rotate_angle=180)
        
        # prepare VASP input files
        for i in toposcarlist:
            prepare_for_VASP(i,name=name)
    
    with open(savepath+'runall.sh','w') as f:
        temp=[]
        for name in namelist:
            temp.append('cd {}_AA && qsub run.sh && cd ..\n'.format(name))
            temp.append('cd {}_AB && qsub run.sh && cd ..\n'.format(name))
        f.writelines(temp)
    
    with open(savepath+'AAAtemp/'+'runall.sh','w') as f:
        temp=[]
        for name in namelist:
            temp.append('cd {}_AAp1 && qsub run.sh && cd ..\n'.format(name))
            temp.append('cd {}_AAp2 && qsub run.sh && cd ..\n'.format(name))
            temp.append('cd {}_ABp1 && qsub run.sh && cd ..\n'.format(name))
            temp.append('cd {}_ABp2 && qsub run.sh && cd ..\n'.format(name))
        f.writelines(temp)


        


