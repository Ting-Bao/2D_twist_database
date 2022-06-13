import argparse, math, numpy as np
from pymatgen.core.structure import Structure
from string import digits, ascii_letters


'''
parser = argparse.ArgumentParser()  # make bilayer from 2D monolayer unit cell
parser.add_argument('--filename', default='POSCAR', type=str,
                    help='Specify monolayer lattice strcture file name')
parser.add_argument('--interlayer', default=3.3, type=float,
                    help='Bilayer interlayer spacing (Ang)')
parser.add_argument('--vacuum', default=12, type=float,
                    help='Specify vacuum layer thickness (Ang)')
parser.add_argument('--ashift', default=0, type=float,
                    help='Top layer shift along lattice vector a (direct)')
parser.add_argument('--bshift', default=0, type=float,
                    help='Top layer shift along lattice vector b (direct)')
parser.add_argument('--rotate_angle', default=180, type=int,
                    help='Top layer rotate angle arround z-axis (degree)')
parser.add_argument('--rotate_anchor', nargs='+', default=[0.5, 0.5, 0], type=float,
                    help='Top layer rotate anchor in direct coordinate of unit cell')
args = parser.parse_args()
'''

def sort_atom_z(str_tmp):
    # str_tmp = pymatgen structure read from POSCAR (vacuum must be along z)
    str_tmp.sort()  # sort coords by element electronegativity
    element_list,element_num = find_element_info(str_tmp)
    frac_coords = str_tmp.frac_coords.copy()
    pn = 0  # number of atoms for previous element
    for i in element_num:
        # sort frac_coords by z-axis values, a[np.argsort(a[:,-1])]
        frac_coords[pn:i+pn] = frac_coords[pn:i+pn][np.argsort(frac_coords[pn:i+pn][:,-1])]
        pn += i
    str_new = Structure(str_tmp.lattice.matrix,
                        str_tmp.species,
                        frac_coords,
                        coords_are_cartesian=False,
                        to_unit_cell=True)
    return str_new

def find_element_info(str_tmp):
    # str_tmp = pymatgen structure read from POSCAR
    # x = file path of POSCAR/*.cif (str)
    # str_tmp = Structure.from_file(x)
    str_formula = str_tmp.formula.split()   # chemical formula
    str_element_name = list(map(lambda str_formula: str_formula.translate(str.maketrans('','',digits)), str_formula))  # store element name
    str_element_num = list(map(lambda str_formula: int(str_formula.translate(str.maketrans('','',ascii_letters))), str_formula))  # store number of each element
    # example: Mn1 Bi2 Te4 = ['Mn','Bi','Te'] & [1,2,4]
    return str_element_name, str_element_num


def make_bilayer(fromposcar,toposcar,interlayer=3.3,vacuum=15,ashift=0,bshift=0,rotate_angle=0,rotate_anchor=[0.5,0,5,0.0]):
    str_tmp = Structure.from_file(fromposcar)  # read in 2D monolayer structure from POSCAR

    # set interlayer spacing and give coords of layer 2
    z_max_layer1 = max(str_tmp.cart_coords[:,2])
    z_min_layer1 = min(str_tmp.cart_coords[:,2])
    str_cart_coords_layer2 = str_tmp.cart_coords.copy()
    str_cart_coords_layer2[:,2] += (interlayer + z_max_layer1 - z_min_layer1)
    str_cart_coords_bilayer = np.vstack((str_tmp.cart_coords,str_cart_coords_layer2))

    # set new lattice constant along z
    str_lattice = str_tmp.lattice.matrix.copy()
    str_lattice[2,2] = math.ceil(vacuum + interlayer + 2*(z_max_layer1 - z_min_layer1))

    str_bilayer = Structure(str_lattice,
                            str_tmp.species*2,
                            str_cart_coords_bilayer,
                            coords_are_cartesian=True,
                            to_unit_cell=True)

    # rotation of 2nd layer
    idx_layer2_start, idx_layer2_end = int(len(str_bilayer)/2), len(str_bilayer)-1  # start <= x < end
    idx_atom_layer2 = np.linspace(idx_layer2_start, idx_layer2_end, idx_layer2_start, dtype=int)  # get second-layer index
    anchor = rotate_anchor[0]*str_bilayer.lattice.matrix[0,:] + rotate_anchor[1]*str_bilayer.lattice.matrix[1,:]
    str_bilayer.rotate_sites(list(idx_atom_layer2),  # list of rotate-site index
                            rotate_angle/180 * np.pi,  # rotation angle in radian
                            np.array([0, 0, 1]),  # rotation axis vector
                            anchor,   # rotation axis anchor
                            to_unit_cell=True)

    # set layer 2 shift along lattice vector a & b
    str_frac_coords_bilayer = str_bilayer.frac_coords.copy()
    str_frac_coords_bilayer[len(str_tmp):len(str_tmp)*2, 0] += ashift
    str_frac_coords_bilayer[len(str_tmp):len(str_tmp)*2, 1] += bshift
    str_bilayer_shift = Structure(str_bilayer.lattice.matrix,
                                str_bilayer.species,
                                str_frac_coords_bilayer,
                                coords_are_cartesian=False,
                                to_unit_cell=True)

    # make sure the atomic coords are arranged properly
    # for each element, first half atoms are in layer 1, while the second in layer 2
    str_bilayer_shift_sort = sort_atom_z(str_bilayer_shift)

    # str_bilayer_shift.sort()
    # fc = str_bilayer_shift.frac_coords.copy()
    # # fc[0,2] = 0.8
    # element_list,element_num = find_element_info(str_bilayer_shift)
    # pn = 0  # previous_atom_num
    # for i in element_num:
    #     fc[pn:i+pn] = fc[pn:i+pn][np.argsort(fc[pn:i+pn][:,-1])]
    #     pn += i
    # str_bilayer_shift_sort = Structure(str_bilayer_shift.lattice.matrix,
    #                                    str_bilayer_shift.species,
    #                                    fc,
    #                                    coords_are_cartesian=False,
    #                                    to_unit_cell=True)

    # output lattice structure to files
    str_bilayer_shift_sort.to(fmt='poscar', filename=toposcar)
