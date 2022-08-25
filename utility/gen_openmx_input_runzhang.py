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

def set_kp(str_tmp):
    # automatically set kpoints (gamma-centered) in SCF calc
    # vacuum layer must be along z
    lattice = str_tmp.lattice.matrix
    lattice_len_a = np.sqrt(lattice[0].dot(lattice[0]))
    lattice_len_b = np.sqrt(lattice[1].dot(lattice[1]))
    lattice_len_c = np.sqrt(lattice[2].dot(lattice[2]))
    kp_a = int(np.ceil(60 / lattice_len_a))
    kp_b = int(np.ceil(60 / lattice_len_b))
    kp_c = int(np.ceil(60 / lattice_len_c))
    return kp_a, kp_b, kp_c

def vec_angle(x1, x2):
    # calculate angle between two vectors
    # x1 & x2 are two vectors (array)
    angle_radian = np.arccos(x1.dot(x2)/np.sqrt(x1.dot(x1))/np.sqrt(x2.dot(x2)))
    angle_degree = angle_radian/np.pi*180
    return angle_radian, angle_degree

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

def openmx_input_norm(struct, basis_accu, soc, dftu, dftuval, magmom, band_calc):
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
    band_path_60 = fr"""
Band.dispersion                 on
Band.Nkpath                     3
<Band.kpath
30 0.500 0.000 0.000  0.000 0.000 0.000  M \Gamma
30 0.000 0.000 0.000  0.666666667 0.333333333 0.000  \Gamma K
30 0.666666667 0.333333333 0.000  0.500 0.000 0.000  K M
Band.kpath>
"""
    band_path_120 = fr"""
Band.dispersion                 on
Band.Nkpath                     3
<Band.kpath
30 0.500 0.000 0.000  0.000 0.000 0.000  M \Gamma
30 0.000 0.000 0.000  0.333333333 0.333333333 0.000  \Gamma K
30 0.333333333 0.333333333 0.000  0.500 0.000 0.000  K M
Band.kpath>
"""
    band_path_90 = fr"""
Band.dispersion                 on
Band.Nkpath                     4
<Band.kpath
30 0.000 0.000 0.000  0.000 0.500 0.000  \Gamma Y
30 0.000 0.500 0.000  0.500 0.500 0.000  Y M
30 0.500 0.500 0.000  0.500 0.000 0.000  M X
30 0.500 0.000 0.000  0.000 0.000 0.000  X \Gamma
Band.kpath>
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

    angle_radian, angle_degree = vec_angle(lattice[0], lattice[1])
    if round(angle_degree) == 60:  # judge lattice angle
        band_key = band_calc * band_path_60
    elif round(angle_degree) == 120:
        band_key = band_calc * band_path_120
    elif round(angle_degree) == 90:
        band_key = band_calc * band_path_90
    else:
        band_key = ''
        print('Lattice symmetry does not match 2D materials')

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
{band_key}
### END ###
"""
    return in_file


pot_path = '/home/xurz/POT/DFT_DATA19'
basis_info_path = '/home/xurz/POT/opmx_basis.txt'

dftuval = ['Mn', 'd', '4.0']
magmom = [5.0] + [-5.0] + [0.0]*12  # example MnBi2Te4 bilayer unit cell

unit_struct, lattice, unit_element_name, unit_element_num, unit_atom_total_num, species, frac_coords, cart_coords = read_struct_2d('POSCAR')
unit_struct = Structure.from_file('POSCAR')
unit_struct.to('poscar', 'POSCAR_crystal')
# openmx_input_norm(struct, basis_accu, soc, dftu, dftuval, magmom, band_calc)
openmx_input = openmx_input_norm(unit_struct, 1, True, True, dftuval, magmom, True)
with open('openmx_in.dat', 'w') as f:
    f.write(openmx_input)
