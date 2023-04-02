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
from pymatgen.core.structure import Structure

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
    # print(len(model))

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
        try:  # in case no all_data.json available
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
            # shutil.copytree(folder_path,temppath)
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

        # print(folder_path)
    print("{} Models found and selected!", format(len(models)))


def infer_and_eval(source, target):
    '''
    get the prepare file to run on the w001 server
    '''
    models = [i for i in os.listdir(source) if os.path.isdir(source + i)]
    to_infer_list = []  # 记录需要比较的cased，eg， 'Al2Cl2/2-1'
    cmd = []  # 记录需要在w001上执行的脚本
    cmd_project = []  # 记录需要在w001上执行的project脚本
    count = 0
    for i in models:
        temppath = target + i

        # check if there corresponding openmx result, else skip
        excludelist = ['P4-276f0a298324']
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

            cmd.append("\necho start processing:{}".format(i + ' ' + twist_case))
            '''
            1. deeph-inference
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
            # print(i)
            cmd.append(
                '\n# cd {} && deeph-inference --config inference.ini && sleep 3 \n'
                .format(temp_workpath))
            '''
            2. e3 eval.py
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
                'cd {} && python {} eval.ini -n 8 >> eval_log && sleep 3 \n'.
                format(temp_workpath + 'predict_e3nn/', e3eval_path))
            # 多进程并行时候，-n 8 在有些case下不work，不知道为何
            '''
            3. get mae
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
                'cd {} && python mae_heat_map.py>>info_mae && sleep 3 \n'.
                format(temp_workpath + 'mae'))
            '''
            4. band_structure and infer_band_pardiso.py 
            '''
            # Use the high symmertric line in k space !!same as that in openmx calcualtion!!
            # Use the Feimi Energy in openmx calculation

            # make dir
            DFT_path = opmx_path_local + '../twisted_band_opmx/{}_{}/'.format(
                i, twist_case)
            if not os.path.exists(twistpath_local + '/band'):
                os.mkdir(twistpath_local + '/band')

            # prepare the pardiso file
            fermi = grep_fermi_opmxout(DFT_path + 'openmx.out')  # Hartree
            kpath = [
                i.replace('\n', '",\n"')
                for i in grep_kpath_from_opmx(DFT_path + 'openmx_in.dat')
            ]
            num_k = sum([int(i.split()[0]) for i in kpath]) + 1
            kpath = '["' + "".join(kpath)[:-4] + '"]'
            soc = grep_json_key(jsonpath=source + i + '/src/dataset_info.json',
                                key='spinful')
            num_band = 140 if twist_case == '2-1' else 280
            template = 'template/analyze/infer_band_pardiso_template.py'
            content = [fermi, kpath, soc, num_k, num_band]
            pardiso = from_template(template=template, content=content)
            with open(twistpath_local + '/band/' + 'infer_band_pardiso.py',
                      'w',
                      encoding='utf-8') as fp:
                fp.writelines(pardiso)

            # get the soft link ln -s
            for file in ['hamiltonians_pred.h5', 'overlaps.h5', 'element.dat', 'lat.dat',
                         'orbital_types.dat', 'rlat.dat', 'R_list.dat', 'site_positions.dat']:
                linkfile = twistpath_local + '/band/{}'.format(file)
                try:
                    os.remove(os.path.abspath(linkfile))
                except:
                    pass
                os.symlink('../predict_e3nn/{}'.format(file), linkfile)
            shutil.copy('template/analyze/plot_band_scatter_each_k.py',
                        twistpath_local + '/band/')

            cmd.append(
                'cd {} && python infer_band_pardiso.py>>info_band && python plot_band_scatter_each_k.py >> info_band&& sleep 3 \n'
                .format(temp_workpath + 'band'))
            '''
            5. pband, band_projection 
            '''
            # make dir
            DFT_path = opmx_path_local + '../twisted_band_opmx/{}_{}/'.format(
                i, twist_case)
            if not os.path.exists(twistpath_local + '/pband'):
                os.mkdir(twistpath_local + '/pband')

            # prepare the pardiso file
            template = 'template/analyze/infer_pband_template.py'
            content = [fermi, kpath, soc, num_k,
                       num_band]  # same as the band part
            pardiso = from_template(template=template, content=content)
            with open(twistpath_local + '/pband/' + 'infer_pband.py',
                      'w',
                      encoding='utf-8') as fp:
                fp.writelines(pardiso)

            template = 'template/analyze/post_orb.sh'
            content = [num_k]
            temp = from_template(template=template, content=content)
            with open(twistpath_local + '/pband/' + 'post_orb.sh',
                      'w',
                      encoding='utf-8') as fp:
                fp.writelines(temp)

            # get the soft link ln -s
            for file in ['hamiltonians_pred.h5', 'overlaps.h5', 'element.dat', 'lat.dat',
                         'orbital_types.dat', 'rlat.dat', 'R_list.dat', 'site_positions.dat']:
                linkfile = twistpath_local + '/pband/{}'.format(file)
                try:
                    os.remove(os.path.abspath(linkfile))
                except:
                    pass
                os.symlink('../predict_e3nn/{}'.format(file), linkfile)

            cmd.append(
                'cd {} && python infer_pband.py>>info_pband && bash post_orb.sh &&sleep 3 \n'
                .format(temp_workpath + 'pband'))
            '''
            6. pdos
            '''
            # make dir
            if not os.path.exists(twistpath_local + '/pdos'):
                os.mkdir(twistpath_local + '/pdos')
            if not os.path.exists(twistpath_local + '/pdos/pdos'):
                os.mkdir(twistpath_local + '/pdos/pdos')

            # get the soft link ln -s
            for file in ['hamiltonians_pred.h5', 'overlaps.h5', 'element.dat', 'lat.dat',
                         'orbital_types.dat', 'rlat.dat', 'R_list.dat', 'site_positions.dat']:
                linkfile = twistpath_local + '/pdos/{}'.format(file)
                try:
                    os.remove(os.path.abspath(linkfile))
                except:
                    pass
                os.symlink('../predict_e3nn/{}'.format(file), linkfile)
            # info.json
            linkfile = twistpath_local + '/pdos/info.json'
            try:
                os.remove(os.path.abspath(linkfile))
            except:
                pass
                os.symlink('../pband/info.json', linkfile)

            # decide kmesh
            str_raw = Structure.from_file(DFT_path + 'POSCAR')
            kmesha = int(120 // str_raw.lattice.a + 1)
            kmeshb = int(120 // str_raw.lattice.b + 1)
            # prepare the dos_config.json
            template = 'template/analyze/dos_config_template.json'
            content = [float(fermi) * 27.21138602, kmesha, kmeshb,
                       num_band]  # same as the band part
            temp = from_template(template=template, content=content)
            with open(twistpath_local + '/pdos/' + 'dos_config.json',
                      'w',
                      encoding='utf-8') as fp:
                fp.writelines(temp)

            f1 = "julia16 /home/baot/bin/sparse_calc_pardis_pdos.jl -i .  -o ./pdos --config dos_config.json"
            f2 = "/home/lihe/anaconda3/envs/3_9/bin/python3 /home/baot/bin/post_pdos_orb_type.py -i1 . -i2 pdos -o pdos"
            f3 = "/home/lihe/anaconda3/envs/3_9/bin/python3 /home/baot/bin/smearing_pdos_orb_type.py -i1 . -i2 pdos -o . --config dos_config.json"
            f4 = "/home/lihe/anaconda3/envs/3_9/bin/python3 /home/baot/bin/plot_pdos_orb_type.py"
            cmd.append("cd {} && ".format(temp_workpath + 'pdos') + f1 + ' \n')
            cmd.append("cd {} && ".format(temp_workpath + 'pdos') + f2 + ' \n')
            cmd.append("cd {} && ".format(temp_workpath + 'pdos') + f3 + ' \n')
            cmd.append("cd {} && ".format(temp_workpath + 'pdos') + f4 + ' \n')

            count += 1
            head = fr'''#!/bin/bash
#PBS -N eval
#PBS -l nodes=1:ppn=64
#PBS -l walltime=960:00:00
#PBS -q cmtmd
source ~/.bashrc
#conda activate 39
'''
            if count % 15 == 0:
                with open(infer_path + 'runall{}.sh'.format(count//15), 'w', encoding='utf-8') as f:
                    f.writelines([head]+cmd)
                cmd = []
    print(" Infer and eval file prepared for {} twisted material(s)! ".format(
        len(to_infer_list)))


if __name__ == '__main__':
    # select_model(source=results_collection, target = model_set)
    # don't run if the model is already refined
    infer_and_eval(source=model_set + 'e3nn_batch/', target=infer_path)
