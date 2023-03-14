import numpy as np 
import os
import argparse
import json
#import matplotlib.pyplot as plt 

parser = argparse.ArgumentParser(description='projected band')
parser.add_argument(
        '-i1','--input_dir_deeph', type=str, default='./',
        help='path of deeph formatted *.dat and *.h5 files'
        )
parser.add_argument(
        '-i2','--input_dir_egval', type=str, default='./',
        help='path of output of julia script of eigenvector, including pband.dat and egval.dat'
        )
parser.add_argument(
        '-o','--output_dir', type=str, default='./',
        help='path of output'
        )
args = parser.parse_args()
input_path = args.input_dir_deeph
input_path2 = args.input_dir_egval
output_path = args.output_dir

# try to build correspondance from global orbital index to (ia, element, n, l, m, spin)
with open("{}/info.json".format(input_path),'r') as info_f:
    info_j = json.load(info_f)
    E_F = float(info_j["fermi_level"])
    isspinful = bool(info_j["isspinful"])

element = np.loadtxt("{}/element.dat".format(input_path))
orb_info = np.empty(0)
with open("{}/orbital_types.dat".format(input_path), 'r') as orb_type_f:
    for ia in range(len(element)):
        #print(ia)
        line = orb_type_f.readline()
        line = line.strip().split()
        for iz in range(len(line)):
            if iz == 0 :
                t_iz = 1
            else:
                if line[iz] == line[iz-1]:
                    t_iz += 1
                else:
                    t_iz = 1
            il = int(line[iz])
            for im in range(2*il+1):
                orb_info = np.append(orb_info, ia+1) # atom index, 1-based
                orb_info = np.append(orb_info, int(element[ia])) # element index
                orb_info = np.append(orb_info, t_iz) # n, 1-based
                orb_info = np.append(orb_info, il) # l, 0-based
                orb_info = np.append(orb_info, im) # m, 0-based
                orb_info = np.append(orb_info, 0) # spin, 0 for spin down
orb_info = np.reshape(orb_info,(int(len(orb_info)/6),6))
np.savetxt("{}/orb_info.dat".format(output_path), orb_info, fmt="%d")
if isspinful:
    t_orb_info = np.empty((len(orb_info)*2,6))
    io = 0
    ia = 1
    t_io = 0
    buff_atom = np.empty(0)
    while io < len(orb_info):
        if ia != orb_info[io,0]:
            buff_atom = np.reshape(buff_atom,(int(len(buff_atom)/6),6))
            for iao in range(len(buff_atom)):
                t_orb_info[t_io, :] = buff_atom[iao, :]
                t_io += 1
            for iao in range(len(buff_atom)):
                t_orb_info[t_io, :] = buff_atom[iao, :]
                t_orb_info[t_io, -1] = 1 # spin up
                t_io += 1
            buff_atom = np.empty(0)
            ia += 1
        buff_atom = np.append(buff_atom,orb_info[io,:])
        if io == len(orb_info)-1:
            buff_atom = np.reshape(buff_atom,(int(len(buff_atom)/6),6))
            for iao in range(len(buff_atom)):
                t_orb_info[t_io, :] = buff_atom[iao, :]
                t_io += 1
            for iao in range(len(buff_atom)):
                t_orb_info[t_io, :] = buff_atom[iao, :]
                t_orb_info[t_io, -1] = 1 # spin up
                t_io += 1
            buff_atom = np.empty(0)
        io += 1

        # if .jld follows what in Hop order ....
        # t_orb_info[:len(orb_info),:] = orb_info.copy()
        # t_orb_info[len(orb_info):,:] = orb_info.copy()
        # t_orb_info[len(orb_info):,-1] = 1
    orb_info = t_orb_info.copy()
    np.savetxt("{}/orb_info.dat".format(output_path), orb_info, fmt="%d")

periodic_table = {0: 'n', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U'}
angular_mom = {0: 's', 1: 'p', 2: 'd', 3: 'f'}

egval = np.loadtxt("{}/egval.dat".format(input_path2))
pband = np.loadtxt("{}/pband.dat".format(input_path2))
# pband shape: (norbit, nband)
norbit = len(orb_info)
nband = len(pband[0])
with open("{}/reduced_pband_orb_type.dat".format(output_path),'w') as out_f:
    for iband in range(nband):
        reduced_pband = {}
        for iorbit in range(norbit):
            #print(iorbit)
            this_key = "{} {}".format(periodic_table[orb_info[iorbit, 1]], angular_mom[orb_info[iorbit, 3]])
            if this_key not in reduced_pband.keys():
                #print(this_key)
                reduced_pband[this_key] = 0
            reduced_pband[this_key] += pband[iorbit,iband]
        #print(egval[iband],reduced_pband)
        out_f.write("{:>12.6f} ".format(egval[iband]-E_F))
        for key in reduced_pband.keys():
            out_f.write("      {:>6s}  {:>10.6f}".format(key, reduced_pband[key]))
        # out_f.write("          S {:>10.6f}".format(reduced_pband['S s']+reduced_pband['S p']+reduced_pband['S d']))
        # out_f.write("          Bi p {:>10.6f}".format(reduced_pband['Bi p']))
        out_f.write("\n")
        
