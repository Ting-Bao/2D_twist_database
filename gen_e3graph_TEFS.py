# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.10
#
# This file is to generate DeepH-E3 graph file and prepare files for Tencent TEFS platform
# This file is supposed to 
#
#    !!!!!!!!!!!!!!!!!!!!!!
#    !!! run on B401GPU !!! 
#    !!!!!!!!!!!!!!!!!!!!!!
#
# if having IO permission problem, try `sudo python thisfile.py` or `sudo chown baot:cmt ./*`



# suggest to use absolute path here, following is path on B401GPU
pythonpath='python'
pfolder='/home/baot/TS216NAS/Public/baot/dataset/processed/' #processed folder
gfolder='/home/baot/TS216NAS/Public/baot/dataset/e3graph/' #where to put the graph
e3process='/home/baot/bin/DeepH-E3-221004_modified/deephe3-preprocess.py'
dataset_note='_3x3_576_xyz_0.1'



import os
import shutil
import json
from utility.utils import *

def read_info_isspinful(file):
    with open(file,'r') as f:
        temp=json.load(f)
        #print(temp["isspinful"])
        #print(type(temp["isspinful"]))
    if temp["isspinful"]:
        return 'soc'
    else:
        return 'nosoc'

materials = os.listdir(pfolder)
#print(materials)
cmd_bash=[]
mat2gen=[]
for i in materials:
    #print(os.path.join(gfolder,i,'graph'))
    if not os.path.exists(os.path.join(gfolder,i,'graph')):
        if not os.path.exists(os.path.join(gfolder,i)):
            os.makedirs(os.path.join(gfolder,i))
        mat2gen.append(i)
        shutil.copy('template/TEFS/train.sh',os.path.join(gfolder,i))
        shutil.copy('template/TEFS/traine3-gpu_spot.json',os.path.join(gfolder,i))
        try:
            soctag=read_info_isspinful(pfolder+i+'/processed/0_0/info.json')
        except:
            soctag=read_info_isspinful(pfolder+i+'/processed/0_0_0/info.json')
        datasetname=i.replace('-', '_')+'_'+soctag+dataset_note
        #write gen_e3graph.ini
        gene3ini=from_template(template='template/DeepH_config/gen_e3graph.ini',\
            content=[os.path.join(pfolder,i,'processed'),os.path.join(gfolder,i,'graph'),datasetname])
        with open(os.path.join(gfolder,i,'gen_e3graph.ini'),'w',encoding='utf-8') as f:
            f.writelines(gene3ini)
        #write train.ini
        trainini=from_template(template='template/TEFS/e3train.ini_template',\
            content=[datasetname])
        with open(os.path.join(gfolder,i,'train.ini'),'w',encoding='utf-8') as f:
            f.writelines(trainini)
        cmd_bash.append('{} {} {} -n 8 \n'.format(pythonpath,e3process,os.path.join(gfolder,i,'gen_e3graph.ini')))



with open(gfolder+'cmd.sh','w',encoding='utf-8') as f:
    f.writelines(cmd_bash)

print('{} materials to gen graph:\n'.format(len(mat2gen)))
for i in mat2gen:
    print(i,end='\n')

if input('The config for training is updated, continue to generate graph? Y/N \n') in ['Y','y']:
    os.system('cd {} && chmod +x cmd.sh && ./cmd.sh'.format(gfolder))
    print('finished!')
else:
    print('exit!')