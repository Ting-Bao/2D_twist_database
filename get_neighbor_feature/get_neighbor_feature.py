import os
import json
import numpy as np
from numpy.linalg import eig
from pymatgen.core.structure import Structure
from argparse import ArgumentParser
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick
import pandas as pd
from scipy.interpolate import interp1d
from functools import partial

def get_parser():
    parser = ArgumentParser(description='')
    parser.add_argument(
        '-f',
        '--filename',
        type=str,
        default=
        'POSCAR',
        help='')
    parser.add_argument('-r', '--rcut', type=float, default=11, help='')
    parser.add_argument('-s', '--sample', type=int, default=50, help='')
    args = parser.parse_args()
    return parser

def log10(x):
    if x>0:
        return np.log10(x)
    else:
        return 0

def get_data(parser,path):
    if os.path.exists(path+'/nb.json'):
        print('Load nb info from existed json file.')
        with open (path+'/nb.json','r',encoding='utf-8') as fp:
            return json.load(fp)
    data = {}
    tw_case = os.listdir(path)
    for i in tw_case:
        if not 'POSCAR' in i:
            tw_case.remove(i)
    data['config'] = get_dataset(parser, path=path + '/config')
    for case in tw_case:
        data[case[7:]] = get_poscar(parser, path=path + '/' + case)
    with open(path+'/nb.json','w',encoding='utf-8') as fp:
        json.dump(data,fp,indent=2)
    return data

def get_poscar(parser, path):
    args = parser.parse_args()
    # get neighboring bond distance of prediction target
    structure_target = Structure.from_file(path)
    nb_target = structure_target.get_all_neighbors(
        r=args.rcut)  # find neighbors
    nb_bond_dist = []  # bond distance of all neighbors
    for atom_each in nb_target:
        for nb_single_idx in atom_each:
            nb_bond_dist.append(nb_single_idx[1])
    return nb_bond_dist


def get_dataset(parser, path):
    # get neighboring bond distance of training set
    args = parser.parse_args()
    train_set_path = path
    train_set_list = os.listdir(train_set_path)
    # compare_set = []
    nb_set_bond_dist_tot = []
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
    return nb_set_bond_dist_tot

def plotfig(parser, data, name):
    args = parser.parse_args()
    colorlist=['red','green','blue']

    # new figure by baot
    fig, ax = plt.subplots(figsize=(9,5))
    #ax2 = ax.twinx()
    ax2 = ax
    # ax2.spines['right'].set_visible(False)
    twist=list(data.keys())
    twist.remove('config')
    twist = sorted(twist,key= lambda x: int(x.split('-')[0]),reverse=False)
    for i in range(len(twist)):
        item = twist[i]
        nb_bond_dist = data[item]
        set_hist, _ = np.histogram(nb_bond_dist,
                               bins=np.linspace(1.0,
                                                args.rcut + 0.2,
                                                args.sample + 1,
                                                endpoint=True))
        # normalized
        # set_hist = set_hist / len(nb_bond_dist) *100
        set_hist=[log10(i) for i in set_hist]
        width = (args.rcut + 0.2 - 1.0) / args.sample
        X = np.linspace(1.0, args.rcut + 0.2, args.sample + 1, endpoint=True)
        Y = np.append(set_hist, 0)
        interp_func = interp1d(X, Y, kind='quadratic')
        num_points = 3 * args.sample
        X_new = np.linspace(1.0, args.rcut + 0.2, num_points,
                            endpoint=False) + width / (2 * 3)
        Y_new = interp_func(X_new)
        ax2.plot(X_new, Y_new, linewidth=1.2, color=colorlist[i])
    ax2.set_ylim(ymin=0)
    # 将坐标轴的base number设置为一位。1是指科学计数法时的位数

    nb_set_bond_dist_tot = data['config']
    set_hist, _ = np.histogram(nb_set_bond_dist_tot,
                               bins=np.linspace(1.0,
                                                args.rcut + 0.2,
                                                args.sample + 1,
                                                endpoint=True))
    # normalized
    # set_hist = set_hist / len(nb_set_bond_dist_tot)*100
    set_hist=[log10(i) for i in set_hist]
    width = (args.rcut + 0.2 - 1.0) / args.sample
    X = np.linspace(1.0, args.rcut + 0.2, args.sample + 1, endpoint=True)
    Y = np.append(set_hist, 0)
    interp_func = interp1d(X, Y, kind='quadratic')
    num_points = 3 * args.sample
    X_new = np.linspace(1.0, args.rcut + 0.2, num_points,
                        endpoint=False) + width / (2 * 3)
    Y_new = interp_func(X_new)

    ax.plot(X_new, Y_new, color='grey', label='dataset')
    tempy=[-0.001] * len(X_new)
    ax.fill_between(X_new, Y_new, tempy, where=(Y_new > tempy), facecolor='grey', alpha=0.5) # set the transparency of the area
    # fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
    # xticks = mtick.FormatStrFormatter(fmt)
    # ax.yaxis.set_major_formatter(xticks)
    def mjrFormatter(x, pos):
        return "$10^{{{0}}}$".format(int(x))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    ax.set_ylabel('Count', color='black')
    ax.set_xlim(xmin=1, xmax=10)
    ax.set_ylim(ymin=0)
    # ax.yaxis.get_major_formatter().set_powerlimits((0, 4))



    # 手动设置legend
    plots = [Line2D([0], [0], color='grey')]+[Line2D([0], [0], color=i) for i in colorlist]
    # labels = ['dataset'] + ['{}'.format(item) for item in twist]
    labels = ['dataset'] + [r'10.87$^\circ$'] + [r'13.17$^\circ$'] + [r'1.08$^\circ$']
    ax.legend(plots, labels, loc='lower right')

    # Set the title and axis labels
    ax.set_title(name.split('-')[0])
    ax.set_xlabel(r'Neighboring atom-atom distance ($\AA$)')

    plt.savefig('get_neighbor_feature/fig/new_statistics.svg')
    plt.savefig('get_neighbor_feature/fig/new_statistics.jpg', dpi=800)
    plt.show()


def main():
    parser = get_parser()
    materials = [
        'C2-a6735a4a3797', 'P4-276f0a298324', 'SnS2-08a9307b286e',
        'Bi2I6-433fccc74b5d', 'MoS2-b3b4685fb6e1', 'PtSe2-d000f0288397'
    ]
    for material in materials[1:2]:
        path = 'get_neighbor_feature/target_POSCAR/' + material
        data = get_data(parser,path)
        
        plotfig(parser, data, name=material)
    print('OK')


if __name__ == '__main__':
    main()