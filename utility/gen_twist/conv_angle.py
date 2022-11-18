# From He Li
from pymatgen.core.structure import Structure
import numpy as np
import argparse


def from_60_to_120(lattice):
    print("convert to 120 from 60 ...")
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    return np.stack([a, b - a, c], axis=0)

def from_120_to_60(lattice):
    print("convert to 60 from 120 ...")
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    return np.stack([a, a + b, c], axis=0)

def conv(source, target, task, fmt='poscar'):
    stru = Structure.from_file(source)

    if task == '60_to_120':
        lattice_new = from_60_to_120(stru.lattice.matrix)
    elif task == '120_to_60':
        lattice_new = from_120_to_60(stru.lattice.matrix)
    elif task == None:
        a_b_angle = stru.lattice.gamma
        if np.allclose(a_b_angle, 120):
            lattice_new = from_120_to_60(stru.lattice.matrix)
        elif np.allclose(a_b_angle, 60):
            lattice_new = from_60_to_120(stru.lattice.matrix)
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError

    stru_new = Structure(
        lattice=lattice_new,
        species=stru.atomic_numbers,
        coords=stru.cart_coords,
        to_unit_cell=True,
        coords_are_cartesian=True
    )
    stru_new.to(fmt, target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_name', type=str, default='POSCAR_120')
    parser.add_argument('--output_name', type=str, default='POSCAR_60')
    parser.add_argument('--task', type=str, default='120_to_60')
    args = parser.parse_args()

    conv(args.input_name, args.output_name, args.task)
