import json, os, deeph

# path naming rule:
# ***/sn_111_novdw_norelax/1-2/H_predict/sn_soc_5x5_17x17x1ln_xyz0.1/all/band_each_k
# put this script in path ***/all/band_each_k
# model_name = 'sn_5x5_30x10rd_xy0.1_z0.2_novdw'  # change when use
train_style = 'all'  # options: all/single/diag3 !!!
# train_spec = 'single'  # training specification !!!
E_Fermi_Hartree = -0.14563531978385  # change when use!!!, grep Chemical openmx.out
lowest_band = -0.21
max_inter = 300
num_band = 80
kp_path = ["15 0.5 0.0 0.0 0.0 0.0 0.0 M Γ","15 0.0 0.0 0.0 0.666666667 0.333333333 0.0 Γ K","15 0.666666667 0.333333333 0.0 0.5 0.0 0.0 K M"]

working_path = os.getcwd()
model_name = working_path.split('/')[-3]
train_spec = working_path.split('/')[-2]
# model_name = os.path.split(os.path.abspath(os.path.join(working_path, '../../')))[-1]
overlap_path = os.path.abspath(os.path.join(working_path, '../../../../only_olp'))
model_path = os.path.abspath(os.path.join(working_path, '../../../../../model', model_name, train_spec))
calc_src_dir = os.path.join(os.path.dirname(deeph.__file__), 'inference', 'sparse_calc.jl')

# write model list
model_path_key = ''
if train_style != 'all':
    for i in os.listdir(model_path):
        model_path_key += f'"{model_path}/{i}",'
    model_path_key = model_path_key[:-1]
else:
    model_path_key = f'"{model_path}"'

# write inference.ini in band prediction working directory
inference_ini = fr"""[basic]
OLP_dir = {overlap_path}
Hop_dir = /home/xurz/bin/Hop.jl
work_dir = {working_path}
structure_file_name = POSCAR_crystal
trained_model_dir = [{model_path_key}]
task = [1, 2, 3, 4, 5]
sparse_calc_config = {working_path}/band_config.json

[interpreter]
julia_interpreter = /home/xurz/bin/julia-1.5.4/bin/julia

[graph]
radius = 7.5
create_from_DFT = True"""
with open(os.path.join(working_path, 'inference.ini'), 'w') as inif:
    inif.write(inference_ini)

# change band_config.json at working dir when use!!!
band_config_no_k_json = {"calc_job": "band",
                         "which_k": -1,
                         "fermi_level": E_Fermi_Hartree*27.21138602,
                         "lowest_band": lowest_band,
                         "max_iter": max_inter,
                         "num_band": num_band,
                         "k_data": kp_path}
with open(os.path.join(working_path, 'band_config.json'), 'w') as jsf:
    json.dump(band_config_no_k_json, jsf)
os.system('deeph-inference --config inference.ini')

# band inference separately for each k point
bash_cmd = ''
for i in range(45):  # index of k points
    # change band_config.json when use!!!
    band_config_json = {"calc_job": "band",
                        "which_k": i+1,
                        "fermi_level": E_Fermi_Hartree*27.21138602,
                        "lowest_band": lowest_band,
                        "max_iter": max_inter,
                        "num_band": num_band,
                        "k_data": kp_path}
    # band_config_json_rectangular = {"calc_job": "band",}
    band_idx = f'{i+1}'
    config_file_path = os.path.join(working_path, 'egval_k', band_idx)
    os.makedirs(config_file_path, exist_ok=True)
    
    with open(os.path.join(config_file_path, 'band_config.json'), 'w') as jsonf:
        json.dump(band_config_json, jsonf)

    bash_cmd += f'julia {calc_src_dir} --input_dir {working_path} --output_dir {config_file_path} --config {config_file_path}/band_config.json &\n'

# submit = fr"""#!/bin/bash
# #PBS -N deeph
# #PBS -l nodes=1:ppn=64
# #PBS -l walltime=96:00:00

# cd ${{PBS_O_WORKDIR}}

# {bash_cmd}
# wait"""
# with open('run.sh','w') as shf:
#     shf.write(submit)

# run sparse calc in parallel style
bash_cmd += 'wait'
os.system(bash_cmd)
