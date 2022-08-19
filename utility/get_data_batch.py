import os, json, argparse
import numpy as np
import math
from pymatgen.core.structure import Structure
from string import digits, ascii_letters
# from xurz_funct import *
###
'''
parser = argparse.ArgumentParser()
parser.add_argument('--basis', default=2, type=int,
                    help='Set basis accuracy, 1=quick/2=standard/3=precise')
parser.add_argument('--magmom', nargs='+', default=[0.0, 0.0, 0.0], type=float,
                    help='Set magnetic moment for each element, [M1,M2,M3]')
parser.add_argument('--soc', default=True, type=bool,
                    help='Enable SOC calculation')
parser.add_argument('--dftu', default=False, type=bool,
                    help='Enable DFT+U calculation')
parser.add_argument('--dftu_value', nargs='+', default=['Mn', 'd', '4.0'], type=str,
                    help='Set DFT + U value')
parser.add_argument('--supercell', nargs='+', default=[3, 3], type=int,
                    help='Make X*Y supercell')
parser.add_argument('--num', default=[24, 24, 1], type=int,
                    help='Number of training structures, [n_a_shift,n_b_shift,n_pert]')
parser.add_argument('--pert', default=[0.1, 0.1, 0.1], type=float,
                    help='Atomic perturbation along x/y/z aixs (Ang)')
args = parser.parse_args()
'''

# main func is at the end of this file
np.random.seed(567)  # set random seed

def cutoff_radius(element_name, accuracy,basisfile):
    r_born_to_ang = 0.52917721  # convert Bohr to Ang
    element_info = find_basis(element_name, accuracy,basisfile=basisfile)
    basis = element_info[3]
    element_and_radius = basis.split('-')[0]
    radius_tmp = filter(lambda ch: ch in '0123456789.', element_and_radius)

    radius_str = ''
    for i in radius_tmp:
        radius_str += i
    radius_float = float(radius_str)

    radius_ang = radius_float * r_born_to_ang * 2
    radius_out = float(math.ceil(2*radius_ang)/2) # this will have a step of 0.5 --Ting BAO

    return radius_out  # in unit of Ang

### output openmx basis settings for specific element
def find_basis(x, b , basisfile):
    # x = element name (str), b = basis accuracy (int)
    # txt file format: Element_name  ZVAL  POT_name  NELECT  Quick  Standard  Precise
    pot_name = x + '_PBE19'
    # with open('opmx_basis.txt', 'r') as opmx_basis:
    with open(basisfile, 'r') as opmx_basis:
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


def get_batch(fromfile='POSCAR',basisfile='/home/xurz/POT/opmx_basis.txt',\
    pot_path = '/home/xurz/POT/DFT_DATA19',out_path = os.path.join('.', 'config'),\
    basis=1,magmom=[0.,0.,0.],soc=True,dftu=False,dftu_value=['Mn', 'd', '4.0'],supercell=[3,3],num=[24,24,1],pert=[0.1,0.1,0.1]):
    """generate opmx input, get the batch for further calculation (finally for the dataset)
    """



    str_unit = Structure.from_file(fromfile)  # read in unit-cell bilayer structure
    unit_formula = str_unit.formula.split()   # chemical formula
    unit_element_name = list(map(lambda unit_formula: unit_formula.translate(str.maketrans('','',digits)), unit_formula))  # store element name
    unit_element_num = list(map(lambda unit_formula: int(unit_formula.translate(str.maketrans('','',ascii_letters))), unit_formula))  # store number of each element
    unit_atom_total_num = sum(list(map(int,unit_element_num)))  # store total number of ions/atoms in unit cell
    # unit_fc = str_unit.frac_coords  # store atomic fractional coordinates
    # print(unit_element_name, unit_element_num, unit_atom_total_num)
    # print(unit_fc)
    # print(unit_fc)

    shift_list = shift_2d_linear(num[0], num[1])
    for shift_index in range(num[0] * num[1]):
        # shift = np.random.rand()-0.5, np.random.rand()-0.5  # interlayer shift use random

        shift = shift_list[shift_index]  # interlayer shift use linear

        unit_fc_tmp = str_unit.frac_coords
        pre_n = 0  # number of atoms/ions for previous element
        for n in unit_element_num:
            # a1,a2 = int(pre_n),int(n/2+pre_n)    # lower & upper bound of layer 1
            b1, b2 = int(n/2+pre_n), int(n+pre_n)  # lower & upper bound of layer 2
            # print(b1,b2)
            unit_fc_tmp[b1:b2,0] += shift[0]
            unit_fc_tmp[b1:b2,1] += shift[1]
            # unit_element_layer1.append(np.reshape(str_unit.species[a1:a2],(n,1)))
            # print(unit_element_layer1)
            # unit_fc_layer2.append(unit_fc[b1:b2])
            # unit_element_layer2.append(str_unit.species[b1:b2])
            pre_n += n
        # print(unit_fc_tmp)
        # unit_fc_layer1 = np.vstack((unit_fc_layer1))
        # unit_element_layer1 = np.vstack((unit_element_layer1))
        # print(unit_element_layer1)
        # unit_fc_layer2 = np.vstack((unit_fc_layer2))
        # unit_element_layer2 = np.vstack((unit_element_layer2))

        # unit_fc_layer2[:,0] += shift_x  # inset perturbation along x
        # unit_fc_layer2[:,1] += shift_y  # inset perturbation along y
        # unit_fc_shift = np.concatenate([unit_fc_layer1, unit_fc_layer1], axis=0)
        # unit_fc_element_list = np.concatenate([unit_element_layer1, unit_element_layer2], axis=0)

        str_shift = Structure(str_unit.lattice,
                            str_unit.species,
                            unit_fc_tmp,
                            coords_are_cartesian=False,
                            to_unit_cell=True)

        str_shift.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, 1]])
        
        for pert_index in range(num[2]):
            super_cart_coords = str_shift.cart_coords

            pert_x = (np.random.rand(len(str_shift), 1) - 0.5) * (pert[0] * 2)
            pert_y = (np.random.rand(len(str_shift), 1) - 0.5) * (pert[1] * 2)
            pert_z = (np.random.rand(len(str_shift), 1) - 0.5) * (pert[2] * 2)
            super_cart_coords += np.concatenate([pert_x, pert_y, pert_z], axis=-1)

            str_shift_pert = Structure(str_shift.lattice,
                                    str_shift.species,
                                    super_cart_coords,
                                    coords_are_cartesian=True,
                                    to_unit_cell=True)

            # generate OpenMX atomic sites
            lattice = str_shift_pert.lattice.matrix
            frac_coords = str_shift_pert.frac_coords
            element = str_shift_pert.species

            # output lattice matrix for openmx
            def_lattice = f'  {lattice[0, 0]:.16f} {lattice[0, 1]:.16f} {lattice[0, 2]:.16f}\n\
    {lattice[1, 0]:.16f} {lattice[1, 1]:.16f} {lattice[1, 2]:.16f}\n\
    {lattice[2, 0]:.16f} {lattice[2, 1]:.16f} {lattice[2, 2]:.16f}'
            # print(def_lattice)

            def_element = ''   # output element info for openmx
            def_dftu_list = '' # output element dftu list for openmx
            list_element = []  # output element list
            list_n_elect = []  # output n_elect list
            for j in str_shift_pert.formula.split():  # for element Mn,Bi,Te
                element_j = j.translate(str.maketrans('','',digits))
                atomic_num, pot_name, n_elect, basis_set = find_basis(element_j, basis, basisfile=basisfile)
                def_element += f'  {element_j}   {basis_set}   {pot_name}\n'  # write element info

                ion_orb_from_basis_txt = list(basis_set.split('-')[1])  # examp [s,3,p,2,d,2,f,1]
                # ion_orb_range = len(ion_orb_from_basis_txt) / 2
                # ion_orb = list(basis_set.split('-')[1])
                ion_orb_name = ion_orb_from_basis_txt[::2]  # terms with odd index
                ion_orb_num = list(map(int,ion_orb_from_basis_txt[1::2]))  # terms with even index
                dftu_str = ''
                for m in ion_orb_name:  # for orbital s,p,d,f
                    for n in range(ion_orb_num[ion_orb_name.index(m)]):  # for orb_occup num
                        if dftu == True and m == dftu_value[1] and n == 0 and element_j == dftu_value[0]:
                            dftu_str += str(n+1) + m + ' ' + dftu_value[2] + ' '
                        else:
                            dftu_str += str(n+1) + m + ' 0.0 '  # write dft+u in order of element info
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
            for i in range(len(str_shift_pert)):
                ion_n_elect = list_n_elect[list_element.index(str(element[i]))]
                ion_magmom = magmom[list_element.index(str(element[i]))]
                s_up = (ion_n_elect + ion_magmom) / 2
                s_dw = (ion_n_elect - ion_magmom) / 2

                coords_set_soc = soc * '0.0  0.0  0.0  0.0  0'
                coords_set_dftu = dftu * 'off'

                frac_coords_str += f'  {i+1}  {str(element[i])}  \
    {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    {s_up}  {s_dw}  {coords_set_soc}  {coords_set_dftu}\n'  # ferromagnetic case

    #             if i < 9:  # anti-ferromagnetic case
    #                 frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_up}  {s_dw}  {coords_set_soc}  {coords_set_dftu}\n'
    #             else:
    #                 frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_dw}  {s_up}  {coords_set_soc}  {coords_set_dftu}\n'

    ###########################################################################################
    #             if args.soc == 0 and args.dftu == 0:
    #                 frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_up}  {s_dw}\n'
                
    #             elif args.soc == 0 and args.dftu == 1:
    #                 frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_up}  {s_dw}  off\n'

    #             elif args.soc == 1 and args.dftu == 1:  ### changed for afm case, be careful
    #                 if i < 9:  # if the case that bottom layer magnetic ion has opposite spin
    #                     frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_dw}  {s_up}  0.0 0.0  0.0 0.0  0  off\n'
    #                 else:
    #                     frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_up}  {s_dw}  0.0 0.0  0.0 0.0  0  off\n'
            
    #             elif args.soc == 1 and args.dftu == 0:
    #                 frac_coords_str += f'  {i+1}  {str(element[i])}  \
    #   {frac_coords[i,0]:.16f}  {frac_coords[i,1]:.16f}  {frac_coords[i,2]:.16f}  \
    #   {s_up}  {s_dw}  0.0 0.0  0.0 0.0  0\n'
    ###########################################################################################

            frac_coords_str = frac_coords_str[:-1]

            ### DFT+U calculation setting
            # U value need to be set in openmx_in.dat
            dftu_calc_setting = fr"""
scf.Hubbard.U                   on
scf.Hubbard.Occupation          dual
<Hubbard.U.values
{def_dftu_list}
Hubbard.U.values>
"""

            index_list = [str(shift_index), str(pert_index)] # folder index, /config/1_1/
            index = '_'.join(index_list)
            folder_path = os.path.join(out_path, index)
            openmx_in_file = os.path.join(folder_path, 'openmx_in.dat')
            os.makedirs(folder_path, exist_ok=True)  # path for saving perturbated structures

            str_shift_pert.to('poscar', os.path.join(folder_path, 'POSCAR_crystal')) # cp structure to POSCAR

            with open(os.path.join(folder_path, 'shift.json'), 'w') as json_f:
                json.dump({'displace_a': shift[0], 'displace_b': shift[1]}, json_f, indent=4) # save perturbation detail


            if soc == True:  # set calc whether SOC
                spin_mode, soc_mode = 'nc', 'on'
            elif soc == False and max(magmom) != 0:
                spin_mode, soc_mode = 'on', 'off'
            else:
                spin_mode, soc_mode = 'off', 'off'

            # if args.dftu == True:  # set calc whether DFT+U
            dftu_key = dftu * dftu_calc_setting
            # else:
            #     dftu_key = ''
            
            kp_a, kp_b, kp_c = set_kp(str_shift_pert) # set SCF gamma-centered kpoints
            kp_key = f'{kp_a}  {kp_b}  1'
            
            in_file = fr"""
System.Name                       openmx
DATA.PATH                         {pot_path}
HS.fileout                        on

Species.Number                    {len(str_shift_pert.formula.split())}
<Definition.of.Atomic.Species
{def_element}
Definition.of.Atomic.Species>
Atoms.Number                      {len(str_shift_pert)}
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

            with open(openmx_in_file, 'w') as f:
                f.write(in_file)
