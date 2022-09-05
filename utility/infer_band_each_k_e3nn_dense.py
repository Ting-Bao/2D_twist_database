import os, json

# use sparse matrix calculation code in deeph
calc_src_dir = '/home/xurz/bin/DeepH-pack/deeph/inference/dense_calc.jl'

# get current working directory
working_path = os.getcwd()

# must-have files in working_path:
# rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, overlaps.h5

# input info for band prediction
E_Fermi_Hartree = -0.159111481720451  # grep Chemical openmx.out !!!
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

band_config_json = {"calc_job": "band",
                    "fermi_level": E_Fermi_Hartree*27.21138602,
                    "k_data": kp_path}
config_path = os.path.join(working_path, 'band_config.json')
with open(config_path, 'w') as jsonf0:
    json.dump(band_config_json, jsonf0)

os.system(f'julia {calc_src_dir} --input_dir {working_path} --output_dir {working_path} --config {config_path}')
