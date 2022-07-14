# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.07
#
# this file generate input file of dataset, the calculation should be done in w001
# 

from utility.utils import *
from utility.get_data_batch import *

namelistpath='data/finalchoice_calculation/part1_prioirty/namelist.json'
#namelistpath='data/finalchoice_calculation/part2_normal/namelist.json'
outfile='data/finalchoice_calculation/part1_prioirty/'
#outfile='data/finalchoice_calculation/part2_normal/'

serverpath='/home/xyz/baot/dataset/' # used on w001, all parts put in this folder together 


def read_soc(preparejson):
    with open(preparejson,'r',encoding='utf-8') as f:
        temp=json.load(f)
    return temp['soc']

def main():
    namelist=read_namelist(namelistpath)
    for name in namelist:
        if not os.path.exists(outfile+name+'/config'):
            soc=read_soc(outfile+name+'/prepare.json')
            get_batch(fromfile=outfile+name+'/'+'relaxed_POSCAR',out_path=outfile+name+'/config/',\
                basisfile='template/OPMX/opmx_basis.txt',soc=soc)
        # 0-575, devide into 0-287, 288-575
        # 0-287
        if True or not os.path.exists(outfile+name+'run_0-287.sh'):
            template='template/OPMX/run_opmx_batch_0-287.sh'
            temp=opmxbatch_w001(template=template,name=name,configpath=serverpath+name)
            with open(outfile+name+'/run_0-287.sh','w',encoding='utf-8') as f:
                f.writelines(temp)
        # 288-575
        if True or not os.path.exists(outfile+name+'run_288-575.sh'):
            template='template/OPMX/run_opmx_batch_288-575.sh'
            temp=opmxbatch_w001(template=template,name=name,configpath=serverpath+name)
            with open(outfile+name+'/run_288-575.sh','w',encoding='utf-8') as f:
                f.writelines(temp)

def gen_run_all():

    print("runall.sh generated!")
            
            
             

if __name__=='__main__':
    main()
    gen_run_all()
