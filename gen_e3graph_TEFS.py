# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.10
#
# This file is to generate DeepH-E3 graph file and prepare files for Tencent TEFS platform
# This file is supposed to 
#   1. run on baot's PC (B403 ubuntu)
#   2. the generated code for graph generating should run on B401GPU 
# if having IO permission problem, try `sudo python thisfile.py` or `sudo chown baot:cmt ./*`


# suggest to use absolute path here
pythonpath='python'
pfolder='/home/baot/TS216NAS/Public/baot/dataset/processed/' #processed folder
gfolder='/home/baot/TS216NAS/Public/baot/dataset/e3graph/' #where to put the graph
prepare='/home/tingbao/Desktop/2D_twist_database/code/template/general/'
e3process='/home/baot/bin/DeepH-E3-221004_modified/deephe3-preprocess.py'


e3g_ini='''
graph_dir =
DFT_data_dir =
processed_data_dir = {}
save_graph_dir = {}
target_data = hamiltonian
dataset_name = Sb2Te3_soc_3x3_256_xyz_0.1
'''.format('112','1241234_sfd')

print(e3g_ini)

import os
import shutil

materials = os.listdir(from_folder)

cmd_bash=[]
mat2gen=[]
for i in materials:
    if not os.path.exists(os.path.join(gfolder,i)):
        os.makedirs(os.path.join(gfolder,i,'gragh'))
        mat2gen.append(i)
        shutil.copy(os.path.join(prepare,'train.sh'))
        cmd_bash.append('{} {} {} -n 8'.format(python_path,e3process,))


with open(gfolder+'cmd.sh','w',encoding='utf-8') as f:
    f.writelines(cmd_bash)

print('materials to gen graph:\n',mat2gen)

if input('continue to generate graph? type Y') in ['Y','y']:
    os.system('cd {} && chmod +x {} && ./{}'.format())