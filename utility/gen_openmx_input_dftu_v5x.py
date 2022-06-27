import argparse
from pymatgen.core.structure import Structure
from string import digits
# from xurz_funct import *
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('--basis', default=1, type=int,
                    help='Set basis accuracy, 1=quick / 2=standard / 3=precise')
parser.add_argument('--bands', default=True, type=bool,
                    help='Enable band structure output')
parser.add_argument('--magmom', nargs='+', default=[0.0, 0.0, 0.0], type=float,
                    help='Set magnetic moment for each element, [M1,M2,M3], can be expanded to elements > 3')
parser.add_argument('--soc', default=False, type=bool,
                    help='Enable SOC calculation')
parser.add_argument('--dftu', default=False, type=bool,
                    help='Enable DFT+U calculation')
parser.add_argument('--dftu_value', default=['Mo','1d','3.0'], type=str,
                    help='Set DFT U value, [element,orbital,U_value]')
args = parser.parse_args()

print(args)
pot_path = '/home/xurz/POT/DFT_DATA19' # potpath should be availble at the running env
#basis_config='/home/xurz/POT/opmx_basis.txt' 
basis_config='/home/tingbao/Desktop/2D_twist_database/code/template/OPMX/opmx_basis.txt' # basis_config can be put locally 

### output openmx basis settings for specific element
def find_basis(x, b):
    # x = element name (str), b = basis accuracy (int)
    # txt file format: Element_name  ZVAL  POT_name  NELECT  Quick  Standard  Precise
    pot_name = x + '_PBE19'
    # with open('opmx_basis.txt', 'r') as opmx_basis:
    with open(basis_config, 'r') as opmx_basis:
        lines = opmx_basis.readlines()
        for line in lines:
            if pot_name in line:
                ZVAL = int(line.split()[1])     # atomic number       (6)
                POT = str(line.split()[2])      # potential name      (C_PBE19)
                NELECT = float(line.split()[3]) # number of electrons (4.0)
                BASIS = str(line.split()[3+b])  # basis setting       (C6.0-s2p2d1)
    return ZVAL, POT, NELECT, BASIS

### automatically set kpoints (gamma-centered) in SCF calc
def set_kp(str_tmp):
    # vacuum layer must be along z
    lattice = str_tmp.lattice.matrix
    lattice_len_a = np.sqrt(lattice[0].dot(lattice[0]))
    lattice_len_b = np.sqrt(lattice[1].dot(lattice[1]))
    lattice_len_c = np.sqrt(lattice[2].dot(lattice[2]))
    kp_a = int(np.ceil(60 / lattice_len_a))
    kp_b = int(np.ceil(60 / lattice_len_b))
    kp_c = int(np.ceil(60 / lattice_len_c))
    return kp_a, kp_b, kp_c

### calculate angle between two vectors
def vec_angle(x1, x2):
    # x1 & x2 are two vectors (array)
    angle_radian = np.arccos(x1.dot(x2)/np.sqrt(x1.dot(x1))/np.sqrt(x2.dot(x2)))
    angle_degree = angle_radian/np.pi*180
    return angle_radian, angle_degree


str_raw = Structure.from_file('POSCAR')  # read in structure
structure = Structure(str_raw.lattice,
                      str_raw.species,
                      str_raw.frac_coords,
                      coords_are_cartesian=False,
                      to_unit_cell=True)  # make sure to_unit_cell
# structure.to('cif', 'crystal.cif')        # cp structure to *.cif
structure.to('poscar', 'POSCAR_crystal')  # cp structure to POSCAR
lattice = structure.lattice.matrix
frac_coords = structure.frac_coords
cart_coords = structure.cart_coords
element = structure.species

# output lattice matrix for openmx
def_lattice = f'  {lattice[0, 0]:.16f} {lattice[0, 1]:.16f} {lattice[0, 2]:.16f}\n\
  {lattice[1, 0]:.16f} {lattice[1, 1]:.16f} {lattice[1, 2]:.16f}\n\
  {lattice[2, 0]:.16f} {lattice[2, 1]:.16f} {lattice[2, 2]:.16f}'
# print(def_lattice)

def_element = ''   # output element info for openmx
def_dftu_list = '' # output element dftu list for openmx
list_element = []  # output element list
list_n_elect = []  # output n_elect list
for j in structure.formula.split():  # for element
    element_j = j.translate(str.maketrans('','',digits))
    atomic_num, pot_name, n_elect, basis_set = find_basis(element_j, args.basis)
    def_element += f'  {element_j}   {basis_set}   {pot_name}\n'

    ion_orb_from_basis_txt = list(basis_set.split('-')[1])  # examp [s,3,p,2,d,2,f,1]
    # ion_orb_range = len(ion_orb_from_basis_txt) / 2
    # ion_orb = list(basis_set.split('-')[1])
    ion_orb_name = ion_orb_from_basis_txt[::2]  # terms with odd index
    ion_orb_num = list(map(int,ion_orb_from_basis_txt[1::2]))  # terms with even index
    dftu_str = ''
    for m in ion_orb_name:  # for orbital s,p,d,f
        for n in range(ion_orb_num[ion_orb_name.index(m)]):  # for orb_occup num
            if element_j == args.dftu_value[0] and str(n+1)+m == args.dftu_value[1]:
                dftu_str += str(n+1) + m + ' ' + args.dftu_value[2] + ' '  # write dftu setting for each ion and orbital
            else:
                dftu_str += str(n+1) + m + ' 0.0 '
    dftu_str = dftu_str[:-1]
    def_dftu_list += f'  {element_j} {dftu_str}\n'

    list_element.append(element_j)
    list_n_elect.append(n_elect)
def_element = def_element[:-1]  # delete the last empty line
def_dftu_list = def_dftu_list[:-1]
# print(list_element)
# print(list_n_elect)
# print(def_element)
# print(def_dftu_list)

frac_coords_str = ''  # output fractional coords for each ion
for i in range(len(str_raw)):
    ion_n_elect = list_n_elect[list_element.index(str(element[i]))]
    ion_magmom = args.magmom[list_element.index(str(element[i]))]
    s_up = (ion_n_elect + ion_magmom) / 2
    s_dw = (ion_n_elect - ion_magmom) / 2

    coords_set_soc = args.soc * '0.0  0.0  0.0  0.0  0'
    coords_set_dftu = args.dftu * 'off'  # orbital polarization

    frac_coords_str += f'  {i+1}  {str(element[i])}  \
  {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
  {s_up}  {s_dw}  {coords_set_soc}  {coords_set_dftu}\n'

#     if args.soc == 0 and args.dftu == 0:
#         frac_coords_str += f'  {i+1}  {str(element[i])}  \
# {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
# {s_up}  {s_dw}\n'
    
#     if args.soc == 0 and args.dftu == 1:
#         frac_coords_str += f'  {i+1}  {str(element[i])}  \
# {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
# {s_up}  {s_dw}  off\n'

#     if args.soc == 1 and args.dftu == 1:
#         frac_coords_str += f'  {i+1}  {str(element[i])}  \
# {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
# {s_up}  {s_dw}  0.0 0.0  0.0 0.0  0  off\n'

frac_coords_str = frac_coords_str[:-1]
# print(frac_coords_str)

### openmx_in.dat template
### band structure calculation setting
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

### DFT+U calculation setting
# U value need to be set in openmx_in.dat
dftu_calc_setting = fr"""
scf.Hubbard.U                   on
scf.Hubbard.Occupation          dual
<Hubbard.U.values
{def_dftu_list}
Hubbard.U.values>
"""


if args.soc == True:
    spin_mode, soc_mode = 'nc', 'on'
elif args.soc == False and max(args.magmom) != 0:
    spin_mode, soc_mode = 'on', 'off'
else:
    spin_mode, soc_mode = 'off', 'off'

dftu_key = args.dftu * dftu_calc_setting

angle_radian, angle_degree = vec_angle(lattice[0],lattice[1])
if round(angle_degree) == 60:  # judge lattice angle
    band_key = args.bands * band_path_60
elif round(angle_degree) == 120:
    band_key = args.bands * band_path_120
elif round(angle_degree) == 90:
    band_key = args.bands * band_path_90
else:
    band_key = ''
    print('Lattice symmetry does not match 2D materials')

kp_a, kp_b, kp_c = set_kp(structure) # set SCF gamma-centered kpoints
kp_key = f'{max(kp_a, 2)}  {max(kp_b, 2)}  1' # 2D case
# kp_key = f'{max(kp_a, 2)}  {max(kp_b, 2)}  {max(kp_c, 2)}' # 3D case

### Main structure for openmx_in.dat
in_file = fr"""
System.Name                       openmx
DATA.PATH                         {pot_path}
HS.fileout                        on

Species.Number                    {len(structure.formula.split())}
<Definition.of.Atomic.Species
{def_element}
Definition.of.Atomic.Species>
Atoms.Number                      {len(structure)}
Atoms.SpeciesAndCoordinates.Unit  FRAC
<Atoms.SpeciesAndCoordinates
{frac_coords_str}
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit            Ang
<Atoms.UnitVectors
{def_lattice}
Atoms.UnitVectors>

scf.XcType                        GGA-PBE   # LDA/LSDA-CA/LSDA-PW/GGA-PBE
scf.ElectronicTemperature         300.0     # default=300 (K) SIGMA in VASP
scf.energycutoff                  300       # default=150 (Ry = 13.6eV)
scf.maxIter                       2000
scf.EigenvalueSolver              Band      # DC/DC-LNO/Krylov/ON2/Cluster/Band
# scf.Ngrid                       
scf.Kgrid                         {kp_key}
scf.criterion                     4e-08     # (Hartree = 27.21eV)
scf.partialCoreCorrection         on

scf.SpinPolarization              {spin_mode}
scf.SpinOrbit.Coupling            {soc_mode}
# scf.Constraint.NC.Spin
# cf.Constraint.NC.Spin.v

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
# scf.Electric.Field                0.0 0.0 1.0    # E-field (1e9 V/m)
# scf.system.charge                 0.0            # carrier doping (e)

# orderN.HoppingRanges              6.0       # default=5.0 (Ang)
# orderN.NumHoppings                2         # default=2
# orderN.KrylovH.order              320
# orderN.KrylovS.order              2000      # default=4*KrylovH.order
# orderN.Expand.Core                on
# orderN.Recalc.Buffer              off
# orderN.Exact.Inverse.S            on

MD.Type                           Nomd      # Nomd (SCF) / NVT_NH (MD)
# MD.maxIter                        5000
# MD.Opt.criterion                  0.0002  # Hartree/Bohr ~ 0.01 eV/A
# MD.TimeStep                       1       # fs
# <MD.TempControl
# 4
# 1     700.0
# 100   700.0
# 400   400.0
# 700   300.0
# MD.TempControl>
# NH.Mass.HeatBath                  22      # mass*Bohr^2, 20 for Carbon

# Dos.fileout                       on
# Dos.Erange                        -10 10
# Dos.Kgrd                        
{band_key}
## Band analysis and Unfolding ##
# Unfolding.Electronic.Band         __Unfolding.Electronic.Band__
# Unfolding.LowerBound              __Unfolding.LowerBound__
# Unfolding.UpperBound              __Unfolding.UpperBound__
# Unfolding.Nkpoint                 __Unfolding.Nkpoint__
# Unfolding.desired_totalnkpt       __Unfolding.desired_totalnkpt__
# <Unfolding.kpoint
# __Unfolding.kpoint__
# Unfolding.kpoint>

# scf.dftD                          __scf.dftD__
# version.dftD                      __version.dftD__

# DFTD3.damp                        __DFTD3.damp__
# DFTD.Unit                         __DFTD.Unit__
# DFTD.rcut_dftD                    __DFTD.rcut_dftD__
# DFTD.cncut_dftD                   __DFTD.cncut_dftD__
# DFTD.IntDirection                 __DFTD.IntDirection__

### END ###
"""

with open('openmx_in.dat', 'w') as f:
    f.write(in_file)
