"""Twist2D Demo."""
"""run in where the file puts"""
#%%
# +------------+
# | Usage Demo |
# +------------+
from twist2d import *
import os
import pandas as pd

def gen_twist(m,n):
    # Create an object for t2d
    twist_demo = Twist2D()

    # Settings for twist supercell generation
    lattice_type = 'rec'
    # 60 triangular = 'tri'   //   120 hexagonal = 'hex'
    # square = 'squ'   //   rectangular = 'rec'

    # set superlattice index m & n
    m = m
    n = n
    # set interlayer distance

    # if lattice_type == 'rec':
    #--> 1st layer
    super_a1_mult = [m, 0]
    super_a2_mult = [0, n]
    twist_demo.add_layer(super_a1_mult, super_a2_mult, layer_dis=2.973962254363104, scs_x=1.6484998174736949, prim_poscar="POSCAR")
    #--> 2nd layer
    super_a1_mult = [m, -1]
    super_a2_mult = [1, n]
    twist_demo.add_layer(super_a1_mult, super_a2_mult, prim_poscar="POSCAR")

    # Twisting the layers
    #  - start_z: The lowest atom's fractional coordinates in z, default 0.1
    #  - super_a3_z: The length of the c vector in z direction, default 20A.
    twist_demo.twist_layers(start_z=0.3194191294526076, super_a3_z=20)

    # Write results to the file
    twist_demo.write_res_to_poscar(filename='POSCAR_{}-{}.vasp'.format(m,n))

    # (Optional) Calculate the twisted angles of each layer in degree 
    twisted_angles = twist_demo.calc_layers_twist_angles()
    print(m, n, twisted_angles)

if __name__=='__main__':
    data = pd.read_csv('Twist_angle_feasible_new.csv')
    oklist=[]
    for i in range(len(data)):
        if data['EST. atom number'][i] < 5000:
            oklist.append((data['m'][i],data['n'][i]))
    for m,n in oklist:
        gen_twist(m=m,n=n)
