# Ting Bao @ Tsinghua University
# Integrated fron Yang Li and He Li's code
#
# This code is used to generate (m,n) twisted structure
# only used in <a,b>=60 degree case
# for 120 degree case, the poscar will turn to 60 first

from .twist2d import *
from pymatgen.core.structure import Structure
import numpy as np
import shutil
import json
# Create an object for t2d


def from_120_to_60(lattice):
    print("convert to 60 from 120 ...")
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    return np.stack([a, a + b, c], axis=0)


def conv(source, target, fmt='poscar'):
    stru = Structure.from_file(source)
    lattice_new = from_120_to_60(stru.lattice.matrix)
    stru_new = Structure(
        lattice=lattice_new,
        species=stru.atomic_numbers,
        coords=stru.cart_coords,
        to_unit_cell=True,
        coords_are_cartesian=True
    )
    stru_new.to(fmt, target)


def check_60(source, target):
    stru = Structure.from_file(source)
    gamma = stru.lattice.gamma
    if gamma-120 < 1e-5:
        conv(source, target)
    else:
        shutil.copy(source, target)


def find_bilayer_dis(poscarpath):
    stru = Structure.from_file(poscarpath)
    zlist = stru.cart_coords[:, 2]
    zlist.sort()
    mid = round(len(zlist)/2)
    return zlist[mid]-zlist[mid-1]


def get_z_length(poscarpath):
    stru = Structure.from_file(poscarpath)
    return stru.lattice.c

# not validated


def move_anchor(poscarpath, jsonpath=None):
    '''json should have two tags:
    mv_a: anchor position in a, using direct coo
    mv_v: anchor position in b, using direct coo
    '''

    # not tested
    if jsonpath != None:
        with open(jsonpath, 'r', encoding='utf-8') as f:
            tmp = json.load(f)
        mv_a = tmp['mv_a']
        mv_b = tmp['mv_b']
        struc = Structure.from_file(poscarpath)
        new_direct_coords = stru.direct_coords
        for i in range(len(new_direct_coords)):
            new_direct_coords[i][0] += new_direct_coords[i][0]+mv_a
            new_direct_coords[i][1] += new_direct_coords[i][1]+mv_b
        stru_new = Structure(
            lattice=struc.lattice,
            species=stru.atomic_numbers,
            coords=new_direct_coords,
            to_unit_cell=True,
            coords_are_cartesian=False
        )
        newpath = os.path.split(poscarpath)[0]+'/poscar_assigned_anchor'
        stru_new.to(fmt='poscar', filename=newpath)
        return newpath
    else:
        return poscarpath

# TODO


def sperate_poscar(poscarpath):
    '''seperate the homo-bilayer structure into 2 poscars, 
    only apply for bilayer in the middle of the poscar
    '''
    print("apply for bilayer in the middle of the poscar")
    stru = Structure.from_file(poscarpath)
    tmp = stru.as_dict()
    half_spices = stru.species[::2]
    zlist = stru.frac_coords[:, 2]
    mid_z = np.median(zlist)
    coo_up = []
    coo_down = []
    for i in stru.frac_coords:
        if i[2] > mid_z:
            coo_up.append(i)
        else:
            coo_down.append(i)

    poscar_up = Structure(
        lattice=stru.lattice,
        species=half_spices,
        coords=coo_up,
        coords_are_cartesian=False
    )
    poscar_up.to(filename=poscarpath+'_up', fmt='poscar')

    poscar_down = Structure(
        lattice=stru.lattice,
        species=half_spices,
        coords=coo_down,
        coords_are_cartesian=False
    )
    poscar_down.to(filename=poscarpath+'_down', fmt='poscar')
    print("bilayer poscar sperated")


def get_twist_struc(m, n, poscar1, poscar2, layer_dis, tofile="POSCAR.T2D.vasp", start_z=0.3, super_a3_z=25.0):
    twist_demo = Twist2D()
    # --> 1st layer
    super_a1_mult = [m, n]
    super_a2_mult = [-n, m+n]
    twist_demo.add_layer(super_a1_mult, super_a2_mult,
                         layer_dis=layer_dis, prim_poscar=poscar1)  # change!!!
    # --> 2nd layer
    super_a1_mult = [n, m]
    super_a2_mult = [-m, n+m]
    twist_demo.add_layer(super_a1_mult, super_a2_mult,
                         prim_poscar=poscar2)  # change!!!

    # Fill the cell with the layers
    # start_z用分数坐标 change!!!
    twist_demo.twist_layers(start_z=0.3, super_a3_z=25)

    # (Optional) Calculate the twisted angles of each layer in degree
    twisted_angles = twist_demo.calc_layers_twist_angles()
    print(twisted_angles)

    # Write results to the file
    twist_demo.write_res_to_poscar(filename=tofile)
    stru = Structure.from_file(tofile)
    atom_num=len(stru.species)
    return twisted_angles,atom_num


def gen_twist(m, n, fromfile='POSCAR', tofile='twistPOSCAR'):
    
    print('gen twisted structure \n from {} to {}'.format(fromfile,tofile))
    
    # STEP 1:check the unitcell has gamma==60 degree first
    check_60(source=fromfile, target=fromfile+'_60')
    file60 = fromfile+'_60'  # this is the path of a make-sure 60 degree poscar file

    # STEP 2:get the layer_distance and z length
    layer_dis = find_bilayer_dis(file60)
    z_lenth = get_z_length(file60)  # return the z_axis length

    # STEP 3: move the anchor
    file60 = move_anchor(poscarpath=file60, jsonpath=None)

    # STEP 4: sperate the bilayer poscar, the atoms should be around the middle of the poscar
    sperate_poscar(poscarpath=file60)

    # STEP 5: gen twisted structure
    twisted_angles, atom_num = get_twist_struc(m=m, n=n, poscar1=file60+'_up', poscar2=file60+'_down', layer_dis=layer_dis,
                    tofile=tofile, super_a3_z=z_lenth)
    
    return twisted_angles, atom_num


if __name__ == '__main__':
    gen_twist(m=3, n=2, fromfile='demo/relaxed_POSCAR',
              tofile='demo/POSCAR_hhh')
    # gen_twist(m=1,n=2)
