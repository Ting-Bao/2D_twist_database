import csv
import numpy as np
# set lattice constant for rectangular lattice
lattice_a = 3.3127999306000002
lattice_b = 4.5033998489000000

# set range for m index
m_range = range(2, 80)
# e.g. (m,0) (0,n) for layer_1, (m,-1) (1,n) for layer 2

# set threshold for strained angle
threshold_ang = 0.6

# initialize variables
twist_ang = 0
ang_list  = []
#ang_str   = ''

for m in m_range:
    twist_vector_a = np.array([m,1])
    twist_a = np.array([twist_vector_a[0]*lattice_a,twist_vector_a[1]*lattice_b])
    len_a = np.sqrt(twist_a[0]**2 + twist_a[1]**2)
    
    for n in range(1, m):
        twist_vector_b = np.array([-1,n])
        twist_b = np.array([twist_vector_b[0]*lattice_a,twist_vector_b[1]*lattice_b])
        len_b = np.sqrt(twist_b[0]**2 + twist_b[1]**2)
        
        strain_ang = np.arccos((twist_a[0]*twist_b[0] + twist_a[1]*twist_b[1]) \
                               /(len_a*len_b)) /np.pi * 180
        
        twist_ang = np.arctan(lattice_b/m/lattice_a) / np.pi * 180
        
        n_atom_est = 4*m*n*2
        
        ang_list.append([m, n, strain_ang, twist_ang, n_atom_est])

ang_list_feasible = []
for i in range(len(ang_list)):
    ang_diff = abs(ang_list[i][2] - 90)
    if ang_diff <= threshold_ang:
        ang_list_feasible.append(ang_list[i])

# headers = ['m','n','strain angle','twist angle','EST. atom number']
# with open('Twist_angle_feasible.csv','w') as f:
#     f_csv = csv.writer(f)
#     f_csv.writerow(headers)
#     f_csv.writerow(ang_list_feasible)

# for i in range(1, np.shape(ang_list)[0]+1):
#     ang_diff = abs(ang_list[i][2] - 90)
#     if ang_list[i][0] == 
        
#        ang_str += f'{m}    {n}    {strain_ang}    {twist_ang}'
#        if i != len(structure_shift_perturb) - 1:
#            ang_str += '\n'
#
#with open('ang_list.dat', 'w') as f:
#    f.write(ang_str)
