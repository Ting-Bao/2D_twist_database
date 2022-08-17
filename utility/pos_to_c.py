import numpy as np
import os

# get all structure file names in current directory
def target_files(str):  # aften use POSCAR here
    working_path = os.getcwd()
    filename_list = []
    for root, dirs, files in os.walk(working_path):  # standard usage for os.walk
        break
    for f in files:
        if str in f:  # get all filenames containing POSCAR
            filename_list.append(f)
    return filename_list

# read in structure from POSCAR file and convert to spglib C input
def struc2c(filename):
    with open(filename, 'r') as f:
        f_lines = f.readlines()
        lattice = np.zeros([3, 3])
        sites = []
        type_atom = []

        for i in range(len(f_lines)):
            if i >= 2 and i <= 4:
                lattice[i-2, :] = list(map(float, f_lines[i].split()))  # spglib input
                
            elif i == 6:
                num_each_specie = list(map(int, f_lines[i].split()))
                num_atom = sum(num_each_specie)  # spglib input
                for j in range(len(num_each_specie)):
                    type_atom += [j+1]*num_each_specie[j]  # spglib input

            elif i >= 8 and i < (8 + num_atom):
                sites_tmp = list(map(float, f_lines[i].split()[0:3]))
                sites.append(sites_tmp)
        
        type_atom = np.array(type_atom) # spglib input
        sites = np.array(sites)  # spglib input
        # four major inputs for spglib: lattice, num_atom, type_atom, sites
        # all in np.array format
    return lattice, type_atom, num_atom, sites

# transfer structure info into string input to spglib C program
def info_spglib_c(lattice, type_atom, num_atom, sites):
    lattice_str, type_atom_str, num_atom_str, sites_str = '', '', '', ''
    
    for i in range(3):
        for j in range(3):
            lattice_str += f'{lattice[j, i]} '  # fixed length of 9
    
    num_atom_str = f'{num_atom} '  # just an int number
    
    for k in range(num_atom):
            type_atom_str += f'{type_atom[k]} '  # fixed length of num_atom
    
    for m in range(num_atom):
        for n in range(3):
            sites_str += f'{sites[m, n]} '  # fixed length of 3*num_atom

    info = lattice_str + num_atom_str + type_atom_str + sites_str
    return info

# main procedure
filename_list = target_files('POSCAR')  # get all filenames containing POSCAR
print(filename_list)
for i in filename_list:
    lattice, type_atom, num_atom, sites = struc2c(i)
    info = info_spglib_c(lattice, type_atom, num_atom, sites)
    print(i)
    print(info)
    in_file = f'{info}'
    file_name = i + '_info.txt'
    with open(file_name, 'w') as f:
        f.write(in_file)

"spglib_2d xxx_info.txt >> "