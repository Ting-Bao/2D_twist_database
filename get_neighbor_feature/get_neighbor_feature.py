import os
import numpy as np
from numpy.linalg import eig
from pymatgen.core.structure import Structure
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from scipy.interpolate import interp1d

def get_parser():
    parser = ArgumentParser(description='')
    parser.add_argument('-f',
                        '--filename',
                        type=str,
                        default='get_neighbor_feature/target_POSCAR/P4-276f0a298324/POSCAR_12-6',
                        help='')
    parser.add_argument('-r', '--rcut', type=float, default=8.0, help='')
    parser.add_argument('-s', '--sample', type=int, default=50, help='')
    args = parser.parse_args()
    return parser

def get_POSCAR(parser,path):
    args = parser.parse_args()
    # get neighboring bond distance of prediction target
    structure_target = Structure.from_file(path)
    nb_target = structure_target.get_all_neighbors(r=args.rcut)  # find neighbors
    nb_bond_dist = []  # bond distance of all neighbors
    for atom_each in nb_target:
        for nb_single_idx in atom_each:
            nb_bond_dist.append(nb_single_idx[1])
    stat_target = np.histogram(
        np.array(nb_bond_dist), args.sample,
        [1.0, args.rcut + 0.2])  # keep x-axis sampling identical
    height_target = stat_target[0]  # only use the histogram height data

def get_dataset(parser,path):
    # get neighboring bond distance of training set
    args = parser.parse_args()
    train_set_path = path
    train_set_list = os.listdir(train_set_path)
    compare_set, nb_set_bond_dist_tot = [], []
    height_set = np.zeros((args.sample))
    for i in train_set_list:
        # read train-set structures
        structure_set = Structure.from_file(
            os.path.join(train_set_path, i, 'POSCAR_crystal'))
        nb_set = structure_set.get_all_neighbors(r=args.rcut)  # find neighbors
        nb_set_bond_dist = []  # bond length of all neighbors
        for atom_set_each in nb_set:
            for nb_set_single_idx in atom_set_each:
                nb_set_bond_dist.append(nb_set_single_idx[1])
                nb_set_bond_dist_tot.append(nb_set_single_idx[1])
        stat_set = np.histogram(np.array(nb_set_bond_dist), args.sample,
                                [1.0, args.rcut + 0.2])
        height_set += stat_set[0]

        height_set = stat_set[0]
        # compare = height_target.dot(height_set) / np.linalg.norm(
            height_target) / np.linalg.norm(height_set)
        # compare_set.append(compare)
    # compare = height_target.dot(height_set)/np.linalg.norm(height_target)/np.linalg.norm(height_set)
    # print(compare)
    # compare_set.append(compare)

    # normalization among all structures from training set
    # compare_normalized = sum(compare_set) / len(compare_set)
    print(compare_normalized)
    return nb_set_bond_dist

def plot(parser,data):
    args = parser.parse_args()
    # height_set_average = height_set_tot / len(compare_set)
    plt.subplot(2, 1, 1)
    plt.hist(nb_set_bond_dist_tot, args.sample, (1.0, args.rcut + 0.2))
    plt.subplot(2, 1, 2)
    plt.hist(nb_bond_dist, args.sample, (1.0, args.rcut + 0.2))
    # plt.hist(height_target, args.sample, (1.0, args.rcut+0.2))
    # plt.show()
    plt.savefig('get_neighbor_feature/fig/sample.jpg', dpi=800)
    plt.close('all')

    # print(len(nb_set_bond_dist_tot)*13*13)


    # new figure by baot
    fig, ax = plt.subplots()

    ax.hist(nb_bond_dist, args.sample, (1.0, args.rcut + 0.2), color='b')
    ax.spines['right'].set_visible(False)
    ax.set_xlim(xmin=1.0, xmax=args.rcut + 0.2)
    ax.set_ylabel('Count', color='b')
    ax.yaxis.get_major_formatter().set_powerlimits((0, 3))
    # 将坐标轴的base number设置为一位。1是指科学计数法时的位数

    ax2 = ax.twinx()
    set_hist, _ = np.histogram(nb_set_bond_dist_tot,
                            bins=np.linspace(1.0,
                                                args.rcut + 0.2,
                                                args.sample + 1,
                                                endpoint=True))
    width = (args.rcut + 0.2 - 1.0) / args.sample

    X = np.linspace(1.0, args.rcut + 0.2, args.sample+1, endpoint=True) 
    Y = np.append(set_hist,0)
    interp_func = interp1d(X, Y, kind='cubic')

    num_points = 3 * args.sample
    X_new = np.linspace(1.0, args.rcut + 0.2, num_points,
                        endpoint=False) + width / (2 * 3)

    Y_new = interp_func(X_new)

    ax2.plot(X_new, Y_new, color='red', label='dataset')
    ax2.set_ylabel('Count', color='r')
    ax2.set_ylim(ymin=0)
    ax2.yaxis.get_major_formatter().set_powerlimits((0, 4))
    # 将坐标轴的base number设置为一位。1是指科学计数法时的位数

    hist_legend = [
        Line2D([0], [0], color='b'),
        Line2D([0], [0], color='b', lw=0, marker='s', markersize=5)
    ]

    plots = [hist_legend[1], Line2D([0], [0], color='red')]
    labels = ['12-6 structure', 'dataset']
    ax.legend(plots, labels, loc='upper left')

    # Set the title and axis labels
    ax.set_title(r'Neighbor atom distance distribution')
    ax.set_xlabel(r'Neighbor atom distance/ $\AA$')

    plt.savefig('get_neighbor_feature/fig/new_statistics.svg', dpi=800)
    plt.savefig('get_neighbor_feature/fig/new_statistics.jpg', dpi=800)
    # Show the figure
    plt.show()

def main():
    parser = get_parser()
    materials =['C2-a6735a4a3797', 'P4-276f0a298324', 'SnS2-08a9307b286e', 'Bi2I6-433fccc74b5d', 'MoS2-b3b4685fb6e1','PtSe2-d000f0288397']
    for material in materials[1:2]:
        data={}
        path = 'get_neighbor_feature/target_POSCAR/' + material
        tw_case = os.listdir(path)
        tw_case.remove('config')
        namelist = [i[7:] for i in tw_case]
        data['config'] = get_dataset(parser,path=path+'/config')
    print('OK')


if __name__=='__main__':
    main()