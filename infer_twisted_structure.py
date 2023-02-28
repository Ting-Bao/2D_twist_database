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
json_alldata = 'data/c2dbdata/jsondata/'
# data/c2dbdata/jsondata/Al2Cl6-b8b35ca6154e.all_data.json

opmx_path_local = './data/twisted_struc_opmx_OK/twisted_band_opmx_processed_with_S/'
infer_path = './data/infer_calculation/'

# @w001
opmx_path = '/home/xurz/temp_baot/twisted_band_opmx/data/twisted_band_opmx_processed_with_S/'
model_path = '/home/xyz/baot/infer_batch/modelset/e3nn_batch/'
work_path = '/home/xyz/baot/infer_batch/infer_calculation/'
e3eval_path = '/home/xyz/bin/DeepH-E3-221004_modified_fixparser/deephe3-eval.py'
pardiso_path = '/home/xyz/bin/infer_band_pardiso_60k.py'


def select_model(source, target):
    '''
    find human-collect models and put into a uniform format
    '''
    models = [i for i in os.listdir(source) if os.path.isdir(source + i)]
    #print(len(model))

    # find materials with mannual adjust json file.
    #   else to auto add the json file
    # finally put the models into e3 and else into model_set path
    for i in models:
        jsonpath = search_file(dirpath=source + i,
                               filename='train_set.json',
                               result_path=[])
        results_path = search_file(dirpath=source + i,
                                   filename='best_model.pkl',
                                   result_path=[])
        try:
            folder_path = os.path.split(results_path[-1])[0]
        except:
            print('!!! Failed case: {}!!!'.format(i))
            continue

        # get properties
        material_name = i.split('-')[0]
        C2DB_id = i.split('-')[1]
        try:  #in case no all_data.json available
            ICSD_id = if_ICSD(json_alldata + i + '.all_data.json',
                              returnid=True)
        except:
            ISCD_id = 'null'
        if ICSD_id == False:
            ICSD_id = 'null'
        SOC = grep_json_key(folder_path + '/src/dataset_info.json', 'spinful')

        print("Start processing {}".format(i))

        if len(jsonpath) == 0:  # no manual json
            temppath = model_set + 'e3nn_batch/' + i
            #shutil.copytree(folder_path,temppath)
            # 上一行做不到覆盖拷贝
            if not os.path.exists(temppath):
                os.mkdir(temppath)
            os.system("cp -r {} {}".format(folder_path + '/*', temppath + '/'))

            template = 'template/DeepH_config/model_tefs.json'
            jsoninfo = from_template(
                template=template,
                content=[material_name, C2DB_id, ICSD_id, SOC])
            with open(temppath + '/train_set.json', 'w',
                      encoding='utf-8') as f:
                f.writelines(jsoninfo)

        # 手动处理下deeph的case
        elif len(jsonpath) > 0:  # has manual json
            nn_type = grep_json_key(jsonpath[0], 'nn_type')
            if nn_type != 'E3':
                print(i)
                raise AttributeError

            e3nn_ver = grep_json_key(jsonpath[0], 'e3nn_ver')
            if e3nn_ver == 'wrong_key':
                e3nn_ver = '0.4.4'
            if e3nn_ver == '0.4.4':
                temppath = model_set + 'e3nn_batch/' + i
            else:
                temppath = model_set + 'old_e3nn/' + i
                print('old:', i)

            if not os.path.exists(temppath):
                os.mkdir(temppath)
            os.system("cp -r {} {}".format(folder_path + '/*', temppath + '/'))
            add_key_to_json(jsonpath=temppath + '/train_set.json',
                            content={"e3nn_ver": e3nn_ver})

        #print(folder_path)
    print("{} Models found and selected!", format(len(models)))


def infer_and_eval(source, target):
    '''
    get the prepare file to run on the w001 server
    '''
    models = [i for i in os.listdir(source) if os.path.isdir(source + i)]
    to_infer_list = []  #记录需要比较的cased，eg， 'Al2Cl2/2-1'
    cmd = []  # 记录需要在w001上执行的脚本

    for i in models:
        temppath = target + i

        # check if there corresponding openmx result, else skip
        excludelist = ['P4-']
        if not os.path.exists(opmx_path_local +
                              '{}_2-1'.format(i)) and i not in excludelist:
            print('skip ', i)
            continue

        if not os.path.exists(temppath):
            os.mkdir(temppath)

        for twist_case in ['2-1', '3-2', '7-4', '9-5']:
            if not os.path.exists(opmx_path_local +
                                  '{}_{}'.format(i, twist_case)):
                continue
                # 有一些3-2因为原子数太多而被排除,所以不需要比较计算, 多数需要排除P4专有的7-4,9-5
            else:
                twistpath_local = os.path.join(temppath, twist_case)
                if not os.path.exists(twistpath_local):
                    os.mkdir(twistpath_local)
                temp_workpath = work_path + i + '/' + twist_case + '/'
                to_infer_list.append('{}/{}'.format(i, twist_case))

            cmd.append("\n echo start processing:{}".format(i + ' ' +
                                                            twist_case))
            '''
            deeph-inference
            '''
            # content 对应inference.ini中的work_dir, OLP_dir, trained_model_dir
            content = [
                temp_workpath + 'predict_e3nn',
                opmx_path + i + '_' + twist_case, model_path + i
            ]
            template = 'template/DeepH_config/inference.ini'
            jsoninfo = from_template(template=template, content=content)
            with open(twistpath_local + '/inference.ini',
                      'w',
                      encoding='utf-8') as f:
                f.writelines(jsoninfo)
            #print(i)
            cmd.append(
                '\n# cd {} && deeph-inference --config inference.ini && sleep 10 \n'
                .format(temp_workpath))
            '''
            e3 eval.py
            '''
            content = [
                model_path + i, temp_workpath, temp_workpath + 'predict_e3nn',
                temp_workpath + 'predict_e3nn'
            ]
            template = 'template/DeepH_config/eval.ini'
            jsoninfo = from_template(template=template, content=content)
            if not os.path.exists(twistpath_local + '/predict_e3nn'):
                os.mkdir(twistpath_local + '/predict_e3nn')
            with open(twistpath_local + '/predict_e3nn/eval.ini',
                      'w',
                      encoding='utf-8') as f:
                f.writelines(jsoninfo)

            cmd.append('cp -r {} {}\n'.format(
                opmx_path + i + '_' + twist_case + '/*',
                temp_workpath + 'predict_e3nn/'))
            cmd.append(
                'cd {} && python {} eval.ini -n 8 >> eval_log && sleep 10 \n'.
                format(temp_workpath + 'predict_e3nn/', e3eval_path))
            # 多进程并行时候，-n 8 在有些case下不work，不知道为何
            '''
            get mae
            '''
            content = [
                temp_workpath + 'predict_e3nn/hamiltonians_pred.h5',
                temp_workpath + 'predict_e3nn', temp_workpath + 'mae/'
            ]
            template = 'template/analyze/mae_heat_map.py'
            jsoninfo = from_template(template=template, content=content)
            if not os.path.exists(twistpath_local + '/mae'):
                os.mkdir(twistpath_local + '/mae')
            with open(twistpath_local + '/mae/mae_heat_map.py',
                      'w',
                      encoding='utf-8') as f:
                f.writelines(jsoninfo)

            cmd.append(
                'cd {} && python mae_heat_map.py>>info_mae && sleep 10 \n'.
                format(temp_workpath + 'mae', e3eval_path))
            '''
            infer_band_pardiso.py 
            '''
            '''
            band_projection 
            '''
            # tangzc didn't finish the soc part

            # TODO
            #cmd.append('cd {}/predict_e3nn && mv ./predict_e3nn/hamiltonians_pred.h5 ./ && python {}\n'.format(temp_workpath,pardiso_path))

    with open(infer_path + 'runall.sh', 'w', encoding='utf-8') as f:
        f.writelines(cmd)
    print(" Infer and eval file prepared for {} twisted material(s)! ".format(
        len(to_infer_list)))


if __name__ == '__main__':
    #select_model(source=results_collection, target = model_set)
    infer_and_eval(source=model_set + '/e3nn_batch/', target=infer_path)
