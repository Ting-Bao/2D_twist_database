# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
#
# this file generate input file of dataset, the calculation should be done in w001 

from utility.utils import *
from utility.get_data_batch import *

namelistpath='data/finalchoice_calculation/part1_prioirty/namelist.json'
outfile='data/finalchoice_calculation/part1_prioirty/'
genrunsh='data/finalchoice_calculation/runall_priority.sh'
mvsh='data/finalchoice_calculation/mv_all_priority.sh'

#namelistpath='data/finalchoice_calculation/part2_normal/namelist.json'
#outfile='data/finalchoice_calculation/part2_normal/'
#genrunsh='data/finalchoice_calculation/runall_normal.sh'
#mvsh='data/finalchoice_calculation/mv_all_normal.sh'

serverpath = '/home/xyz/baot/dataset/' # used on w001, all parts put in this folder together 
storepath = '/home/xyz/nas_disk/baot/dataset/'
local_basisfile='template/OPMX/opmx_basis.txt'

def read_soc(preparejson):
    with open(preparejson,'r',encoding='utf-8') as f:
        temp=json.load(f)
    return temp['soc']

def get_ele_name(poscarpath):
    with open(poscarpath,'r',encoding='utf-8') as f:
        temp=f.readlines()
    return temp[5].split()

def main():
    namelist=read_namelist(namelistpath)
    for name in namelist:
        if not os.path.exists(outfile+name+'/config'):
            soc=read_soc(outfile+name+'/prepare.json')
            get_batch(fromfile=outfile+name+'/'+'relaxed_POSCAR',out_path=outfile+name+'/config/',\
                basisfile='template/OPMX/opmx_basis.txt',soc=soc)
        # 0-575
        if True or not os.path.exists(outfile+name+'run_opmx.sh'):
            template='template/OPMX/run_opmx_batch.sh'
            temp=from_template(template=template,content=[name,serverpath+name,storepath])
            with open(outfile+name+'/run_opmx.sh','w',encoding='utf-8') as f:
                f.writelines(temp)

def gen_preprocessini():
    namelist=read_namelist(namelistpath)
    for name in namelist:
        template='template/DeepH_config/preprocess.ini'
        radius = max([cutoff_radius(i,1,basisfile=local_basisfile) for i in get_ele_name(outfile+name+'/relaxed_POSCAR')])
        temp=from_template(template=template,content=[serverpath+name,str(radius)])
        with open(outfile+name+'/preprocess.ini','w',encoding='utf-8') as f:
            f.writelines(temp)

def gen_graphini():
    """ use epoches = 0 to gen the graph
    """
    namelist=read_namelist(namelistpath)
    template='template/DeepH_config/gen_graph.ini'
    for name in namelist:
        soctag='soc' if read_soc(outfile+name+'/prepare.json') else 'nosoc'
        temp=from_template(template=template,content=[serverpath+name+'/graph',serverpath+name+'/result',serverpath+name+'/processed',name.replace('-','_')+'_'+soctag+'3x3_24x24x1_xyz0.1'])
        with open(outfile+name+'/gen_graph.ini','w',encoding='utf-8') as f:
            f.writelines(temp)

def gen_mvtonas():
    namelist=read_namelist(namelistpath)
    template='template/general/mv_to_nas.sh'
    for name in namelist:
        temp=from_template(template=template,content=[name,serverpath+name,storepath])
        with open(outfile+name+'/mv_to_nas.sh','w',encoding='utf-8') as f:
            f.writelines(temp)

def gen_run_all():
    namelist=read_namelist(namelistpath)
    cmd=[]
    mv_to_nas=[]
    for name in namelist:
        temp1='cd '+serverpath+name+" && qsub run_opmx.sh\n"
        cmd.append(temp1)
        temp2='cd '+serverpath+name+" && chmod +x mv_to_nas.sh && ./mv_to_nas.sh\n"
        mv_to_nas.append(temp2)
    with open(genrunsh,'w',encoding='utf-8') as f:
        f.writelines(cmd)
    with open(mvsh,'w',encoding='utf-8') as f:
        f.writelines(mv_to_nas)
    print("runall.sh and mv_to_nas.sh generated!")
    
if __name__=='__main__':
    main()
    gen_preprocessini()
    gen_graphini()
    gen_mvtonas()
    gen_run_all()
