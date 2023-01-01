import json
import numpy as np
import h5py
import os
from matplotlib.pylab import plt

rhp = h5py.File("/home/xurz/work/moire-twist/mos2/mos2_ab_vdw12/6-7/predict_e3nn/hamiltonians_pred.h5","r",)
# rhp = h5py.File("/home/xurz/work/moire-twist/snse/7-6/predict_e3nn/snse_aa_3x3_24x24x1ln_xyz0.1_vdw12/e3nn_all_batch1_lr0.002/band_each_k/hamiltonians_pred.h5","r",)

processed_structure = '/home/xurz/work/moire-twist/mos2/mos2_ab_vdw12/6-7/processed' # containing hamiltonians.h5, element.dat, orbital_types.dat, info.json

# attention to save_path
# simplified_output = False

rh = h5py.File(os.path.join(processed_structure, 'hamiltonians.h5'), "r")
element = np.loadtxt(os.path.join(processed_structure, 'element.dat')).astype(int)
orbital_types_dir = os.path.join(processed_structure, 'orbital_types.dat')
with open(os.path.join(processed_structure, 'info.json'), 'r') as f:
    spinful = json.load(f)['isspinful']

orbital_types = []
with open(orbital_types_dir) as f:
    line = f.readline()
    while line:
        orbital_types.append(list(map(int, line.split())))
        line = f.readline()
atom_num_orbitals = [sum(map(lambda x: 2 * x + 1, atom_orbital_types)) * (1 + spinful) for atom_orbital_types in orbital_types]

hopping_keys = []
for atom_i in element:
    for atom_j in element:
        hopping_key = f'{atom_i} {atom_j}'
        if hopping_key not in hopping_keys:
            hopping_keys.append(hopping_key)
    
dif = {}
for hopping_key in hopping_keys:
    atom_i, atom_j = hopping_key.split()
    index_i = np.where(element==int(atom_i))[0][0]
    index_j = np.where(element==int(atom_j))[0][0]
    dif[hopping_key] = np.full((atom_num_orbitals[index_i], atom_num_orbitals[index_j]), 0.0)

dif_num = {}
for key in dif.keys():
    dif_num[key] = 0
    
for key in rhp.keys():
    atom_i = eval(key)[3] - 1
    atom_j = eval(key)[4] - 1
    N_M_str = f'{element[atom_i]} {element[atom_j]}'
    dif[N_M_str] += np.abs(np.array(rh[key]) - np.array(rhp[key]))
    dif_num[N_M_str] += 1

dif_avg = {}
dif_min = {}
dif_max = {}
for key in dif.keys():
    dif[key] /= dif_num[key]
    dif_avg[key] = np.mean(dif[key])
    dif_min[key] = np.amin(dif[key])
    dif_max[key] = np.amax(dif[key])

for k, v in dif.items():
    print('----------')
    print(k)
    print('max:', dif_max[k])
    print('avg:', dif_avg[k])
    print('min:', dif_min[k])
    mae_heatmap = np.zeros([v.shape[0], v.shape[1]])
    tag_retrain = []
    num_retrain = 0
    # if not simplified_output:
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            # print(f'{v[i, j]:8.6f}', end=" ")
            mae_heatmap[i, j] = float(v[i, j])  # runzhang
            # if v[i, j] > 0.005:   # runzhang
                # hop_tag = '{"'+str(k)+'": '+'['+str(i)+', '+str(j)+']}'
                # tag_retrain.append(hop_tag)
                # num_retrain += 1
        # print()

    filename = '-'.join(str(k).split())
    np.savetxt((filename+'.txt'), mae_heatmap*1000)

    plt.figure()
    plt.imshow(mae_heatmap*1000, cmap=plt.cm.Blues, vmin=0, vmax=dif_max[k]*1200)  # runzhang
    plt.colorbar(label='MAE (meV)')
    plt.xlabel('orbital i')
    plt.ylabel('orbital j')
    plt.savefig(filename, dpi=500)

    # print('need retrain:', tag_retrain)  # runzhang
    # print('need retrain #:', num_retrain)


# mode = "mae"  # choices: max / mse / mae
# if mode == "mse":
#     difs = np.full((38, 38), 0.0)
#     num = 0
#     for key in rhp.keys():
#         difs += np.power(np.abs(np.array(rh[key]) - np.array(rhp[key])),2)
#         num += 1
#     dif = difs / num
# elif mode == "max":
#     difs = []
#     for key in rhp.keys():
#         difs.append(np.abs(np.array(rh[key]) - np.array(rhp[key])))
#     difs = np.stack(difs, axis=0)
#     dif = np.amax(difs, axis=0)
# elif mode == "mae":
#     difs = np.full((38, 38), 0.0)
#     num = 0
#     for key in rhp.keys():
#         difs += np.abs(np.array(rh[key]) - np.array(rhp[key]))
#         num += 1
#     dif = difs / num


# for i in range(dif.shape[0]):
#     for j in range(dif.shape[1]):
#         print(dif[i, j], end=" ")
#     print()

rh.close()
rhp.close()
