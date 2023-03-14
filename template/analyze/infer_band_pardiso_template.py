import os, json
# use sparse matrix calculation code in deeph
# calc_src_dir = '/home/xurz/bin/DeepH-pack/deeph/inference/sparse_calc.jl'
calc_src_dir = '/home/baot/bin/sparse_calc_pardis.jl'

os.system("source /home/lihe/intel/oneapi/setvars.sh")

# get current working directory
working_path = os.getcwd()

# must-have files in working_path:
# rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, overlaps.h5

# input info for band prediction
E_Fermi_Hartree = CONTENT1  # grep Chemical openmx.out !!!
lowest_band = -2.05 # no use in pardiso case
max_inter = 300
num_band = CONTENT5

# kp_path = ["15 0.000 0.000 0.000  0.000 0.500 0.000  Γ Y",
#            "15 0.000 0.500 0.000  0.500 0.500 0.000  Y M",
#            "15 0.500 0.500 0.000  0.500 0.000 0.000  M X",
#            "15 0.500 0.000 0.000  0.000 0.000 0.000  X Γ"]  # square or rectangular lattice

# kp_path = ["30 0.5 0.0 0.0 0.0 0.0 0.0 M Γ",
#            "30 0.0 0.0 0.0 0.666666667 0.333333333 0.0 Γ K",
#            "30 0.666666667 0.333333333 0.0 0.5 0.0 0.0 K M"]  # triangular lattice

kp_path = CONTENT2
header = f'source /home/lihe/intel/oneapi/setvars.sh\nexport OMP_NUM_THREADS=64\nexport JULIA_NUM_THREADS=64\nexport LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}:/home/lihe/local/usr/lib64:/home/lihe/lib\n'

# enable for soc
with open('info.json', 'w') as f:
    spinful = {"isspinful": CONTENT3} #False/True
    json.dump(spinful, f)

bash_cmd = ''
#for i in [2]:
for i in range(CONTENT4):  # index of k points, change when use !!!
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
        os.system(header + f'julia16 {calc_src_dir} --input_dir {working_path} --output_dir {working_path} --config {config_path}')
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

        bash_cmd += f'julia16 {calc_src_dir} --input_dir {working_path} --output_dir {config_file_path} --config {config_file_path}/band_config.json &\n'
        # os.system(f'julia {calc_src_dir} --input_dir {working_path} --output_dir {config_file_path} --config {config_file_path}/band_config.json')

        # use the full multi-core processor
        if i%10==0: 
            bash_cmd += '\nwait\n'

bash_cmd += '\nwait\n'
# with open('diag.sh','w',encoding='utf-8') as f:
#     f.write(header + bash_cmd)
os.system(header + bash_cmd)
print("finished!")
