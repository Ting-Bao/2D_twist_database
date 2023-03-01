import os, numpy as np
from numpy.linalg import eig
from pymatgen.core.structure import Structure
from argparse import ArgumentParser
import matplotlib.pyplot as plt

parser = ArgumentParser(description='')
parser.add_argument('-f', '--filename', type=str, default='POSCAR_12-6', help='')
parser.add_argument('-r', '--rcut', type=float, default=7.0, help='')
parser.add_argument('-s', '--sample', type=int, default=500, help='')
args = parser.parse_args()

# get neighboring bond distance of prediction target
structure_target = Structure.from_file(args.filename)
nb_target = structure_target.get_all_neighbors(r = args.rcut) # find neighbors
nb_bond_dist = [] # bond distance of all neighbors
for atom_each in nb_target:
    for nb_single_idx in atom_each:
        nb_bond_dist.append(nb_single_idx[1])
stat_target = np.histogram(np.array(nb_bond_dist), args.sample, [1.8,  args.rcut+0.2]) # keep x-axis sampling identical
height_target = stat_target[0] # only use the histogram height data


# get neighboring bond distance of training set
train_set_path = os.path.join('.', 'config')
train_set_list = os.listdir(train_set_path)
compare_set, nb_set_bond_dist_tot = [], []
height_set = np.zeros((args.sample))
for i in train_set_list:
    # read train-set structures
    structure_set = Structure.from_file(os.path.join(train_set_path, i, 'crystal.cif'))
    nb_set = structure_set.get_all_neighbors(r = args.rcut) # find neighbors
    nb_set_bond_dist = [] # bond length of all neighbors
    for atom_set_each in nb_set:
        for nb_set_single_idx in atom_set_each:
            nb_set_bond_dist.append(nb_set_single_idx[1])
            nb_set_bond_dist_tot.append(nb_set_single_idx[1])
    stat_set = np.histogram(np.array(nb_set_bond_dist), args.sample, [1.8, args.rcut+0.2])
    height_set += stat_set[0]

    height_set = stat_set[0]
    compare = height_target.dot(height_set)/np.linalg.norm(height_target)/np.linalg.norm(height_set)
    compare_set.append(compare)
# compare = height_target.dot(height_set)/np.linalg.norm(height_target)/np.linalg.norm(height_set)
# print(compare)
# compare_set.append(compare)

# normalization among all structures from training set
compare_normalized = sum(compare_set)/len(compare_set)
print(compare_normalized)

# height_set_average = height_set_tot / len(compare_set)
plt.subplot(2,1,1)
plt.hist(nb_set_bond_dist_tot, args.sample, (1.8, args.rcut+0.2))
plt.subplot(2,1,2)
plt.hist(nb_bond_dist, args.sample, (1.8, args.rcut+0.2))
# plt.hist(height_target, args.sample, (1.8, args.rcut+0.2))
# plt.show()
plt.savefig('sample.jpg')

print(len(nb_set_bond_dist_tot)*13*13)