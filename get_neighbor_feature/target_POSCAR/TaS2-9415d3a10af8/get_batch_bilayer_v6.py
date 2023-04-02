import os, json, argparse
import numpy as np
from pymatgen.core.structure import Structure
from string import digits, ascii_letters

### output openmx basis settings for specific element
def find_basis(x, b, p):
    # x = element name (str), b = basis accuracy (int), p = basis info path (txt file)
    # txt file format: Element_name  ZVAL  POT_name  NELECT  Quick  Standard  Precise
    pot_name = x + '_PBE19'
    # with open('opmx_basis.txt', 'r') as opmx_basis:
    with open(p, 'r') as opmx_basis:
        lines = opmx_basis.readlines()
        for line in lines:
            if pot_name in line:
                ZVAL = int(line.split()[1])     # atomic number       (6)
                POT = str(line.split()[2])      # potential name      (C_PBE19)
                NELECT = float(line.split()[3]) # number of electrons (4.0)
                BASIS = str(line.split()[3+b])  # basis setting       (C6.0-s2p2d1)
    return ZVAL, POT, NELECT, BASIS

### make 2d interlayer linear shift coords list
def shift_2d_linear(step_a, step_b):
    # number of shifts along lattice vector a & b
    shift_a, shift_b = np.linspace(0, 1, step_a+1)[0:-1], np.linspace(0, 1, step_b+1)[:-1]
    shift = []
    for i in shift_a:
        for j in shift_b:
            shift.append([i,j])
    return shift

### automatically set kpoints (gamma-centered) in SCF calc
def set_kp(struct):
    # vacuum layer must be along z, read in structure from pymatgen
    lattice = struct.lattice.matrix
    lattice_len_a = np.sqrt(lattice[0].dot(lattice[0]))
    lattice_len_b = np.sqrt(lattice[1].dot(lattice[1]))
    lattice_len_c = np.sqrt(lattice[2].dot(lattice[2]))
    kp_a = int(np.ceil(60 / lattice_len_a))
    kp_b = int(np.ceil(60 / lattice_len_b))
    kp_c = int(np.ceil(60 / lattice_len_c))
    return kp_a, kp_b, kp_c

def read_struct_2d(poscar):
    # read in POSCAR and store required info by DeepH
    struct = Structure.from_file(poscar) # unit cell
    unit_struct = Structure(struct.lattice, struct.species, struct.frac_coords, coords_are_cartesian=False, to_unit_cell=True)
    unit_formula = unit_struct.formula.split()   # get chemical formula, example: Bi2Se3 = Bi2 Se3
    unit_element_name = list(map(lambda unit_formula: unit_formula.translate(str.maketrans('','',digits)), unit_formula)) # ['Bi', 'Se']
    unit_element_num = list(map(lambda unit_formula: int(unit_formula.translate(str.maketrans('','',ascii_letters))), unit_formula)) # 2 3
    unit_atom_total_num = sum(list(map(int, unit_element_num))) # 2+3=5
    lattice = unit_struct.lattice.matrix
    frac_coords = unit_struct.frac_coords
    cart_coords = unit_struct.cart_coords
    species = unit_struct.species # list, ['Bi', 'Bi', 'Se', 'Se', 'Se']
    return unit_struct, lattice, unit_element_name, unit_element_num, unit_atom_total_num, species, frac_coords, cart_coords

### experimental, so far only use for homo-multilayer
def locate_vdw_layer(struct, num_layer):
    # locate atom index within each vdw layer, read in pymatgen structure and number of vdw layers
    # output (2d list): for 2nd atom in 1st layer use vdw_layer_list[0][1]
    cart_coords = struct.cart_coords
    cart_coords_z = cart_coords[:, 2]
    height_max, height_min = np.max(cart_coords_z) + 0.5, np.min(cart_coords_z) - 0.5
    height_all_layer = height_max - height_min
    height = height_all_layer / num_layer
    vdw_layer_list = []
    for i in range(num_layer):
        limit_dn = height_min + i * height
        limit_up = height_min + (i + 1) * height
        layer_idx_1 = np.argwhere(cart_coords_z > limit_dn)
        layer_idx_2 = np.argwhere(cart_coords_z < limit_up)
        layer_idx = list(np.intersect1d(layer_idx_1, layer_idx_2))
        vdw_layer_list.append(layer_idx)
    return vdw_layer_list

# new locate_vdw_layer with list input of num of atoms within each vdw layer (from bottom to top)
def locate_vdw_layer_new(struct, num_atom_each_layer):
    cart_coords = struct.cart_coords
    cart_coords_z = cart_coords[:, 2]
    idx_coords_z = np.argsort(cart_coords_z) # return original idx of elements in sorted list (from small to large)
    vdw_layer_list = []
    pre = 0
    for i in num_atom_each_layer:
        idx_each_vdw_layer = list(idx_coords_z[pre : pre + i])
        vdw_layer_list.append(idx_each_vdw_layer)
        pre += i
    return vdw_layer_list

def make_pert_struct(struct, pert):
    # read in structure from pymatgen
    # pert (list): max magnitude of random atomic perturbation along x/y/z
    cart_coords = struct.cart_coords
    pert_x = (np.random.rand(len(struct), 1) - 0.5) * (pert[0] * 2)
    pert_y = (np.random.rand(len(struct), 1) - 0.5) * (pert[1] * 2)
    pert_z = (np.random.rand(len(struct), 1) - 0.5) * (pert[2] * 2)
    cart_coords += np.concatenate([pert_x, pert_y, pert_z], axis=-1)
    struct_after_pert = Structure(struct.lattice, struct.species, cart_coords, coords_are_cartesian=True, to_unit_cell=True)
    return struct_after_pert

def make_shift_struct(struct, shift, type, num_atom_each_layer, shift_which_layer):
    # read in structure from pymatgen
    # shift (list): shift of frac coords along lattice a/b
    # type (int): interlayer shift, 1=linera, 2=random
    # num_layer (int): num of vdw layers
    # shift_which_layer (int): idx of layer to shift, starting from 0
    if type == 1:
        shift_fc = shift  # linear
    elif type == 2:
        shift_fc = np.random.rand()-0.5, np.random.rand()-0.5  # random
    # vdw_layer_list = locate_vdw_layer(struct, num_layer) # ion idx of each vdw layer
    vdw_layer_list = locate_vdw_layer_new(struct, num_atom_each_layer) # ion idx of each vdw layer, change input from "num_layer" to "num_atom_each_layer"
    unit_fc_tmp = struct.frac_coords
    for n in vdw_layer_list[shift_which_layer]:
        unit_fc_tmp[n, 0] += shift_fc[0]
        unit_fc_tmp[n, 1] += shift_fc[1]
    struct_after_shift = Structure(struct.lattice, struct.species, unit_fc_tmp, coords_are_cartesian=False, to_unit_cell=True)
    return struct_after_shift, shift_fc

def openmx_input(struct, basis_accu, soc, dftu, dftuval, magmom):
    lattice = struct.lattice.matrix
    frac_coords = struct.frac_coords
    element = struct.species
    openmx_lattice = f'  {lattice[0, 0]:.16f} {lattice[0, 1]:.16f} {lattice[0, 2]:.16f}\n\
                         {lattice[1, 0]:.16f} {lattice[1, 1]:.16f} {lattice[1, 2]:.16f}\n\
                         {lattice[2, 0]:.16f} {lattice[2, 1]:.16f} {lattice[2, 2]:.16f}'
    openmx_element = ''   # output element info for openmx
    openmx_dftu_list = '' # output element dftu list for openmx
    list_element = []  # output element list
    list_n_elect = []  # output n_elect list
    chem_formula = struct.formula.split()
    for j in chem_formula: # loop for each element
        element_j = j.translate(str.maketrans('','',digits)) # get element name
        atomic_num, pot_name, n_elect, basis_set = find_basis(element_j, basis_accu, basis_info_path)
        openmx_element += f'  {element_j}   {basis_set}   {pot_name}\n'  # write element info
        
        ion_orb_from_basis_txt = list(basis_set.split('-')[1])  # example [s,3,p,2,d,2,f,1]
        ion_orb_name = ion_orb_from_basis_txt[::2]  # example [s,p,d,f]
        ion_orb_num = list(map(int,ion_orb_from_basis_txt[1::2]))  # example [3,2,2,1]
        dftu_str = ''
        for m in ion_orb_name:  # for orbital s,p,d,f
            for n in range(ion_orb_num[ion_orb_name.index(m)]):  # for orb_occup num
                if dftu == True and m == dftuval[1] and n == 0 and element_j == dftuval[0]:
                    dftu_str += str(n+1) + m + ' ' + dftuval[2] + ' '
                else:
                    dftu_str += str(n+1) + m + ' 0.0 '  # write dft+u in order of element info
        dftu_str = dftu_str[:-1]
        openmx_dftu_list += f'  {element_j} {dftu_str}\n'  # write DFT+U info

        list_element.append(element_j)
        list_n_elect.append(n_elect)

    openmx_frac_coords_list = ''
    for i in range(len(struct)):  # loop for each atom/ion
        ion_n_elect = list_n_elect[list_element.index(str(element[i]))]
        ion_magmom = magmom[i]  # set magmom for each ion
        s_up = (ion_n_elect + ion_magmom) / 2
        s_dw = (ion_n_elect - ion_magmom) / 2
        coords_set_soc = soc * '0.0  0.0  0.0  0.0  0'
        coords_set_dftu = dftu * 'off'
        openmx_frac_coords_list += f'  {i+1}  {str(element[i])}  {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  {s_up}  {s_dw}  {coords_set_soc}  {coords_set_dftu}\n'

    # delete the last empty line
    openmx_element = openmx_element[:-1]
    openmx_dftu_list = openmx_dftu_list[:-1]
    openmx_frac_coords_list = openmx_frac_coords_list[:-1]

    dftu_calc_setting = fr"""
scf.Hubbard.U                   on
scf.Hubbard.Occupation          dual
<Hubbard.U.values
{openmx_dftu_list}
Hubbard.U.values>
"""

    if soc == True:  # set calc whether SOC or magnetic
        spin_mode, soc_mode = 'nc', 'on'
    elif soc == False and max(magmom) > 0:
        spin_mode, soc_mode = 'on', 'off'
    else:
        spin_mode, soc_mode = 'off', 'off'

    dftu_key = dftu * dftu_calc_setting  # set calc whether DFT+U
    kp_a, kp_b, kp_c = set_kp(struct) # set SCF gamma-centered kpoints
    kp_key = f'{max(kp_a, 2)}  {max(kp_b, 2)}  1'

    in_file = fr"""
System.Name                       openmx
DATA.PATH                         {pot_path}
HS.fileout                        on
Species.Number                    {len(chem_formula)}
<Definition.of.Atomic.Species
{openmx_element}
Definition.of.Atomic.Species>
Atoms.Number                      {len(struct)}
Atoms.SpeciesAndCoordinates.Unit  FRAC
<Atoms.SpeciesAndCoordinates
{openmx_frac_coords_list}
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit            Ang
<Atoms.UnitVectors
{openmx_lattice}
Atoms.UnitVectors>
scf.XcType                        GGA-PBE   # LDA/LSDA-CA/LSDA-PW/GGA-PBE
scf.ElectronicTemperature         300.0     # default=300 (K) SIGMA in VASP
scf.energycutoff                  300       # default=150 (Ry = 13.6eV)
scf.maxIter                       2000
scf.EigenvalueSolver              Band      # DC/DC-LNO/Krylov/ON2/Cluster/Band
# scf.Ngrid                       
scf.Kgrid                         {kp_key}
scf.criterion                     4e-06     # (Hartree = 27.21eV)
scf.partialCoreCorrection         on
scf.SpinPolarization              {spin_mode}
scf.SpinOrbit.Coupling            {soc_mode}
# scf.Constraint.NC.Spin
# scf.Constraint.NC.Spin.v
scf.Mixing.Type                   RMM-DIISK
scf.Init.Mixing.Weight            0.3
scf.Mixing.History                30
scf.Mixing.StartPulay             6
scf.Mixing.EveryPulay             1
{dftu_key}
1DFFT.EnergyCutoff                3600
1DFFT.NumGridK                    900
1DFFT.NumGridR                    900
scf.ProExpn.VNA                   off
MD.Type                           Nomd      # Nomd (SCF) / NVT_NH (MD)
### END ###
"""
    return in_file

def make_all_super_shift_pert_struct(poscar, super, shift, type, num_pert, pert, num_atom_each_layer, shift_which_layer):
    # poscar (str): POSCAR file path
    # super (list): period of supercell along lattice vector [3, 3]
    # shift (list): number of shift along lattice vector [16, 16]
    # num_pert: number of random atomic perturbation
    # pert (list): perturbation magnitude along x/y/z axis [0.1, 0.1, 0.1]
    # type (int): type of interlayer shift, 1/2, 1=linear, 2= random

    # read in unit-cell bilayer structure
    unit_struct, lattice, unit_element_name, unit_element_num, unit_atom_total_num, species, frac_coords, cart_coords = read_struct_2d(poscar)
    # make list of interlayer shift along a/b in frac_coords
    shift_fc_list = shift_2d_linear(shift[0], shift[1])
    struct_list, shift_list = [], []
    # loop over interlayer shift
    # for shift_idx in range(shift[0] * shift[1]):
    for shift_idx_1 in range(shift[0]):
        for shift_idx_2 in range(shift[1]):
            idx = shift_idx_2 + shift_idx_1 * shift[1]
            struct_shift, shift_fc = make_shift_struct(unit_struct, shift_fc_list[idx], type, num_atom_each_layer, shift_which_layer) # shift linearly 2nd layer of bilayer
            struct_shift.make_supercell(([[super[0], 0, 0], [0, super[1], 0], [0, 0, 1]]))
            for pert_idx in range(num_pert):
                struct_shift_pert = make_pert_struct(struct_shift, pert)
                struct_list.append(struct_shift_pert)
                shift_list.append(shift_fc)
    return struct_list, shift_list

def batch_file_output(poscar, super, shift, type, num_pert, pert, num_atom_each_layer, shift_which_layer, basis_accu, soc, dftu, dftuval, magmom):
    struct_after_super_shift_pert_list, shift_list = make_all_super_shift_pert_struct(poscar, super, shift, type, num_pert, pert, num_atom_each_layer, shift_which_layer)
    # print(len(struct_after_super_shift_pert_list))
    for shift_idx_1 in range(shift[0]):
        for shift_idx_2 in range(shift[1]):
            for pert_idx in range(num_pert):
                idx_list = [str(shift_idx_1), str(shift_idx_2), str(pert_idx)]
                idx = '_'.join(idx_list)
                folder_path = os.path.join('.', 'config', idx) # folder naming rule, shift1_shift2_pert, example: 5_5_0
                openmx_in_file = os.path.join(folder_path, 'openmx_in.dat')
                os.makedirs(folder_path, exist_ok=True) # path for saving files

                struct_idx = pert_idx + shift_idx_2 * num_pert + shift_idx_1 * shift[1]
                print(struct_idx)
                struct = struct_after_super_shift_pert_list[struct_idx] # get structure for each idx
                struct.to('poscar', os.path.join(folder_path, 'POSCAR_crystal'))

                with open(os.path.join(folder_path, 'shift.json'), 'w') as shift_f:
                    json.dump({'displace_a': shift_list[struct_idx][0], 'displace_b': shift_list[struct_idx][1]}, shift_f, indent=4)

                in_file = openmx_input(struct, basis_accu, soc, dftu, dftuval, magmom)
                with open(openmx_in_file, 'w') as f:
                    f.write(in_file)

pot_path = '/home/xurz/POT/DFT_DATA19'
basis_info_path = '/home/tingbao/Desktop/2D_twist_database/code/get_neighbor_feature/opmx_basis.txt'
np.random.seed(42)  # set random seed

# batch_file_output(poscar, super, shift, type, 
#                   num_pert, pert, num_layer, shift_which_layer, 
#                   basis_accu, soc, dftu, dftuval, magmom) # use locate_vdw_layer
# batch_file_output(poscar, super, shift, type, 
#                   num_pert, pert, num_atom_each_layer, shift_which_layer, 
#                   basis_accu, soc, dftu, dftuval, magmom) # use locate_vdw_layer_new
dftuval = ['Mn', 'd', '4.0']
magmom = [0.0]*9*6  # example MnBi2Te4 3x3 bilayer supercell
# batch_file_output('POSCAR', [3, 3], [16, 16], 1,
#                   1, [0.15, 0.15, 0.15], 2, 1,
#                   1, True, True, dftuval, magmom)
batch_file_output('POSCAR', [3, 3], [24, 24], 1,
                  1, [0.1, 0.1, 0.1], [3, 3], 1,
                  1, False, False, dftuval, magmom)