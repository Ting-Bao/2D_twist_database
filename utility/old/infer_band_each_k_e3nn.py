import os, json

# use sparse matrix calculation code in deeph
calc_src_dir = '/home/xurz/bin/DeepH-pack/deeph/inference/sparse_calc.jl'

# get current working directory
working_path = os.getcwd()

# must-have files in working_path:
# rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, overlaps.h5

# input info for band prediction
E_Fermi_Hartree = -0.159111481720451  # grep Chemical openmx.out !!!
lowest_band = -0.51
max_inter = 300
num_band = 50
# kp_path = ["15 0.000 0.000 0.000  0.000 0.500 0.000  Γ Y",
#            "15 0.000 0.500 0.000  0.500 0.500 0.000  Y M",
#            "15 0.500 0.500 0.000  0.500 0.000 0.000  M X",
#            "15 0.500 0.000 0.000  0.000 0.000 0.000  X Γ"]  # square or rectangular lattice
kp_path = ["15 0.5 0.0 0.0 0.0 0.0 0.0 M Γ",
           "15 0.0 0.0 0.0 0.666666667 0.333333333 0.0 Γ K",
           "15 0.666666667 0.333333333 0.0 0.5 0.0 0.0 K M"]  # triangular lattice

# enable for soc
with open('info.json', 'w') as f:
    spinful = {"isspinful": True}
    json.dump(spinful, f)

bash_cmd = ''
for i in range(45+1):  # index of k points, change when use !!!
    if i == 0:
        band_config_json = {"calc_job": "band",
                            "which_k": -1,
                            "fermi_level": E_Fermi_Hartree*27.21138602,
                            "lowest_band": lowest_band,
                            "max_iter": max_inter,
                            "num_band": num_band,
                            "k_data": kp_path}
        config_path = os.path.join(working_path, 'band_config.json')
        with open(config_path, 'w') as jsonf0:
            json.dump(band_config_json, jsonf0)
        os.system(f'julia {calc_src_dir} --input_dir {working_path} --output_dir {working_path} --config {config_path}')
        print('Finish get jld and Band file')
    else:
        band_config_json = {"calc_job": "band",
                            "which_k": i,
                            "fermi_level": E_Fermi_Hartree*27.21138602,
                            "lowest_band": lowest_band,
                            "max_iter": max_inter,
                            "num_band": num_band,
                            "k_data": kp_path}
        band_idx = f'{i}'
        config_file_path = os.path.join(working_path, 'egval_k', band_idx)
        os.makedirs(config_file_path, exist_ok=True)

        with open(os.path.join(config_file_path, 'band_config.json'), 'w') as jsonf:
            json.dump(band_config_json, jsonf)

        bash_cmd += f'julia {calc_src_dir} --input_dir {working_path} --output_dir {config_file_path} --config {config_file_path}/band_config.json &\n'

bash_cmd += 'wait'
os.system(bash_cmd)

