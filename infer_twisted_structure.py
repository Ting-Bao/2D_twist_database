# author:   Ting Bao @ Tsinghua University
# email:    bao-ting@foxmail.com
# date:     2022.12
#
# this file is condcuted after calc_twisted_structure.py, 
#   which is used for compare the MAE of Hamiltonian matrix element and band structrue
#
# This file should be run on 403ubuntu
# train_set.json will be renamed as info.json

from utility.utils import *
import os
import shutil
import json

results_collection = '/home/tingbao/Desktop/TS216NAS/Public/baot/results_collection/'
model_set = '/home/tingbao/Desktop/TS216NAS/Public/baot/modelset/'
json_alldata='data/c2dbdata/jsondata/'
# data/c2dbdata/jsondata/Al2Cl6-b8b35ca6154e.all_data.json

opmx_path = './data/twisted_band_opmx/'
model_path = model_set
infer_path = './data/infered_results/'

def select_model(source,target):
    '''
    find human-collect models and put into a uniform format
    '''
    models=[ i for i in os.listdir(source) if os.path.isdir(source+i) ]
    #print(len(model))
    
    # find materials with mannual adjust json file. 
    #   else to auto add the json file 
    # finally put the models into e3 and else into model_set path
    for i in models:
        jsonpath=search_file(dirpath=source+i,filename='train_set.json',result_path=[])
        results_path=search_file(dirpath=source+i,filename='best_model.pkl',result_path=[])
        try:
            folder_path=os.path.split(results_path[-1])[0]
        except:
            print('!!! Failed case: {}!!!'.format(i))
            continue
        if len(jsonpath)==0: # no manual json 
            continue
            temppath = model_set+'e3nn_batch/'+i
            #shutil.copytree(folder_path,temppath)
            # 上一行做不到覆盖拷贝
            os.system("cp -r {} {}".format(folder_path,temppath))

            template='template/DeepH_config/model_tefs.json'
            # get properties
            material_name=i.split('-')[0]
            C2DB_id=i.split('-')[1]
            ICSD_id = if_ICSD(json_alldata+i+'.all_data.json',returnid=True)
            if ICSD_id==False:
                ICSD_id='null'
            SOC=grep_json_key(temppath+'/src/dataset_info.json', 'spinful')

            jsoninfo=from_template(template=template,content=[material_name,C2DB_id,ICSD_id,SOC])
            with open(temppath+'/train_set.json','w',encoding='utf-8') as f:
                f.writelines(jsoninfo)
        
        # 手动处理下deeph的case
        elif len(jsonpath)>0:
            nn_type=grep_json_key(jsonpath[0],'nn_type')
            e3nn_ver=grep_json_key(jsonpath[0],'e3nn_ver')
            print(i,e3nn_ver)
            if nn_type!='E3':
                print(i)
                raise AttributeError
        #print(folder_path)



    print("{} Models found!",format(len(models)))




if __name__=='__main__':
    select_model(source=results_collection, target = model_set)

